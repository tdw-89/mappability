# Proposal: Nextflow Pipeline on AWS Batch (Seqera Platform) for Zebrafish Broad-Mark Mappability

## 1. Objective & Context

We want to compute empirical paired-end mappability across **all 28,407 protein-coding genes** in the zebrafish reference genome (`GRCz12tu`, `GCF_049306965.2`) using `UniMap.jl`, tailored to the read lengths and fragment-length distributions of the 87 broad-peak ChIP (`H3K9me2`, `H3K9me3`) and `Input` sequencing runs in `/Users/thomaswolfe/temp_for_analysis/fish/zebrafish/mappability_broad_marks/broad_peak_marks.csv`.

A sequential local bash reference implementation exists at:
- `/Users/thomaswolfe/temp_for_analysis/fish/zebrafish/mappability_broad_marks/run_sim.sh`

Because full-genome protein-coding simulation generates **~2.76 billion read pairs** (~1,000 vCPU-hours of `bowtie2` alignment), running sequentially on an 8-thread workstation would take **~4.5–5 days**. Moving this workflow to a 4-process **Nextflow pipeline executed on AWS Batch via Seqera Platform** reduces end-to-end wall time to **~3.5–4 hours** at an estimated Spot cost of **~$20**.

---

## 2. Input Files & Simulation Parameters

### Local Reference & Assay Metadata Paths
Located under `/Users/thomaswolfe/temp_for_analysis/fish/zebrafish/mappability_broad_marks/`:
- **Reference FASTA**: `ref_links/GCF_049306965.2_GRCz12tu_genomic.fna`
- **FASTA Index**: `ref_links/GCF_049306965.2_GRCz12tu_genomic.fna.fai`
- **Annotation GFF3**: `ref_links/GCF_049306965.2_GRCz12tu_genomic.gff`
- **Bowtie2 Index**: `ref_links/bowtie2/GRCz12tu.{1..4,rev.1,rev.2}.bt2` (prefix: `GRCz12tu`)
- **Assay summary table**: `broad_peak_marks.csv`

### Chosen Grid Parameters
Derived from [simulation_findings.md](./simulation_findings.md) and the quantile distributions in `broad_peak_marks.csv`:

| Parameter | Value | Rationale |
|---|---|---|
| **Target genes** | `28,407` protein-coding genes | GFF3 features with `$3 == "gene"` and `gene_biotype=protein_coding` (excludes `lncRNA`, `pseudogene`, `tRNA`, `rRNA`, etc.). |
| **Flanks** | `--upstream 1500 --downstream 1500` | Default strand-aware flanks; spans 1,251.8 Mb across 28,407 loci, merging into **12,256 simulation regions** (1,060.0 Mb) at $F_{\max}=500$. |
| **Read lengths (`-L`)** | `60,150` | Safe median read lengths for `Input` (`60`) and broad histone marks `H3K9me2`/`H3K9me3` (`150`). |
| **Fragment lengths (`-F`)** | `50,75,100,125,150,175,200,250,300,350,400,450,500` | Full 13-point grid. In `broad_peak_marks.csv`, short Input runs ($q_{50}=86\text{–}97\text{ bp}$) place up to 21.8% of coverage weight on $F=75$, while long H3K9me3/me2 runs ($q_{50}=249\text{–}268\text{ bp}, q_{90}=399\text{–}412\text{ bp}$) place up to 13.3% of coverage weight on $F=400$. All 87 runs have `fragment_length_max = 500`. |
| **Tiling step (`--step`)** | `10` | Halves read count vs `--step 5` with minimal error ($\text{p99}|\Delta T| = 0.024$) and satisfies $\text{step} \le L_{\min}/5 = 12$. |
| **Orientation (`--orientation`)** | `plus` | Matches `both` at $\text{MAPQ} \ge 2$ while halving read count. |
| **Bowtie2 max insert (`-X`)** | `500` | Must match $F_{\max} = 500$. |
| **MAPQ cutoffs (`--cutoffs`)** | `0,2,5,10,20,30` | Default cutoffs covering all `bowtie2` MAPQ plateaus (`{0,1}`, `{2,3,5}`, `{8,10,20,30}`). |

---

## 3. Critical `UniMap.jl` Invariants for Nextflow / AWS S3

The agent implementing the Nextflow pipeline **must** account for the following six invariants in `UniMap.jl`:

1. **Sharding `bowtie2` across batches is mathematically exact**:
   - `bowtie2` aligns each read pair independently using fixed `-X 500` bounds and seeds its tie-breaking PRNG deterministically from `(read_name, seq, qual)` ([simulation_findings.md](./simulation_findings.md), Section 5).
   - Every simulated read encodes its origin in its QNAME (`regionID|chrom|fragStart0|F|orient`; see [tiler.jl](../src/tiler.jl)).
   - Therefore, aligning FASTQ batches on separate Spot instances produces identical alignments to a single monolithic `bowtie2` run.

2. **Never sort or filter BAMs — concatenate shards with `samtools cat`**:
   - [scan_bam!](../src/score.jl) streams `aln/sim_L{L}.bam` in a single pass and requires `R1` and `R2` of each pair to be immediately adjacent (`bowtie2` output order).
   - Regions/batches may appear in any order in the BAM because each pair is looked up by parsing `regionID` from the read name.
   - **Do not** run `samtools sort`, coordinate `samtools merge`, or `--no-unal`. Combine shard BAMs for a given read length `L` using **`samtools cat -o sim_L${L}.bam shard_L${L}_*.bam`**, which concatenates raw BGZF blocks in ~1–2 minutes without recompression and preserves mate adjacency.

3. **Prevent S3 small-file explosion in `SIMULATE` (49,024 FASTQ files $\to$ 80 batch files)**:
   - `unimap simulate` writes one `_R1.fq.gz` and `_R2.fq.gz` per merged region per read length: $12,256 \times 2 \times 2 = 49,024$ files.
   - Staging 49,024 small files to/from an S3 `work/` bucket will severely bottleneck Nextflow.
   - Because concatenating `.fq.gz` files with `cat` yields a valid gzip stream, the `SIMULATE` process should concatenate the per-region FASTQs (in matching `R1`/`R2` order!) into $B \approx 40$ balanced batch pairs per read length (`batch_L{L}_{01..40}_R{1,2}.fq.gz`, ~800 MB per pair) and **`rm -rf sim/fastq`** before the process exits.

4. **`unimap score` does NOT read FASTQ files**:
   - [score_simulation](../src/score.jl) re-enumerates simulated fragments directly from the reference FASTA, verifies counts against `tiling.tsv`, and streams `aln/sim_L{L}.bam`.
   - The `SCORE` process only needs:
     - The metadata files from `SIMULATE`: `simulate.toml`, `loci.tsv`, `regions.tsv`, `tiling.tsv` (~5 MB total)
     - The reference FASTA + `.fai`
     - The concatenated `sim_L60.bam` and `sim_L150.bam`

5. **Pass `-f <fasta>` explicitly to `unimap score`**:
   - `simulate.toml` records the absolute path of the FASTA inside the `SIMULATE` task's container working directory.
   - Because `SCORE` runs in a separate Nextflow task directory, `unimap score` **must** be invoked with `-f "$FASTA"` so it does not look for the stale path recorded in `simulate.toml`.

6. **Tar `mappability.zarr` before S3 staging**:
   - `scores/mappability.zarr` is a Zarr v2 directory with `CHUNK_POSITIONS = 16_384` ([score.jl](../src/score.jl)). Across 1,251.8 Mb of loci and 2 read lengths, it contains **~450,000 chunk files** (~6.5 GB total, already `zstd`-compressed by Blosc).
   - To avoid hundreds of thousands of S3 `PutObject` calls when publishing results, tar the Zarr directory without recompression (`tar -cf scores/mappability.zarr.tar -C scores mappability.zarr && rm -rf scores/mappability.zarr`) inside the `SCORE` task before staging out.

---

## 4. Proposed Nextflow Pipeline Architecture

```
GFF3 + FASTA (.fai)
        │
        ▼
┌──────────────────────┐
│  1. BUILD_LOCI       │  (2 vCPU, 4 GB RAM, ~30s)
│  awk + unimap loci   │
└──────────┬───────────┘
           │ loci.tsv (28,407 genes)
           ▼
┌──────────────────────┐
│ 2. SIMULATE_AND_BATCH│  (32 vCPU, 32 GB RAM, ~15-20m, ~100 GB local scratch)
│ unimap simulate      │  Writes sim metadata + 49,024 region FASTQs, then
│ + cat into 40 batches│  concatenates into 40 batches/L and deletes sim/fastq/
└────┬────────────┬────┘
     │            │
     │ sim_meta/  │ 80 tuples: (L, batch_id, batch_R1.fq.gz, batch_R2.fq.gz)
     │            ▼
     │   ┌──────────────────────┐
     │   │ 3. ALIGN_BATCH (x80) │  (16 vCPU, 16 GB RAM, Spot, ~25-35m each)
     │   │ bowtie2 -X 500 |     │  80 independent tasks (40 batches x 2 read lengths)
     │   │ samtools view -b     │
     │   └──────────┬───────────┘
     │              │ 80 shard BAMs grouped by L
     ▼              ▼
┌──────────────────────┐
│ 4. CAT_AND_SCORE     │  (4 vCPU, 16 GB RAM, ~2.5-3h, ~150 GB local scratch)
│ samtools cat +       │  Concatenates shards into sim_L60.bam & sim_L150.bam,
│ unimap score + tar   │  runs unimap score, tars mappability.zarr
└──────────────────────┘
```

### Process Details

#### Process 1: `BUILD_LOCI`
- **Resources**: `cpus 2`, `memory '4 GB'`, `time '15m'`
- **Inputs**: `path fasta`, `path fai`, `path gff`
- **Outputs**: `path "loci.tsv"`
- **Script**:
  ```bash
  awk -F'\t' '
      $1 !~ /^#/ && $3 == "gene" && $9 ~ /(^|;)gene_biotype=protein_coding(;|$)/ {
          n = split($9, attrs, ";")
          for (i = 1; i <= n; i++) {
              if (attrs[i] ~ /^ID=/) {
                  print substr(attrs[i], 4)
                  break
              }
          }
      }
  ' "$gff" > protein_coding_genes.txt

  julia -t ${task.cpus} --project=/opt/UniMap.jl -m UniMap loci \
      -f "$fasta" \
      -a "$gff" \
      -g protein_coding_genes.txt \
      -o loci.tsv
  ```
  *(Note: Extracting `ID=` rather than `Name=` is required because 23 protein-coding gene symbols occur multiple times in `GRCz12tu` and would otherwise be skipped as ambiguous aliases).*

#### Process 2: `SIMULATE_AND_BATCH`
- **Resources**: `cpus 32`, `memory '32 GB'`, `disk '120 GB'`, `time '2h'`
- **Inputs**: `path fasta`, `path fai`, `path loci_tsv`
- **Params**: `read_lengths = "60,150"`, `fragment_lengths = "50,75,100,125,150,175,200,250,300,350,400,450,500"`, `step = 10`, `orientation = "plus"`, `n_batches = 40`
- **Outputs**:
  - `path "sim_meta/*"`, emit: `meta` (`simulate.toml`, `loci.tsv`, `regions.tsv`, `tiling.tsv`)
  - `path "batches/batch_L*_R1.fq.gz"`, `path "batches/batch_L*_R2.fq.gz"` (paired into channel tuples `(L, batch_id, r1, r2)`)
- **Script**:
  ```bash
  julia -t ${task.cpus} --project=/opt/UniMap.jl -m UniMap simulate \
      -f "$fasta" \
      -l "$loci_tsv" \
      -o sim \
      -L "${params.read_lengths}" \
      -F "${params.fragment_lengths}" \
      --step ${params.step} \
      --orientation ${params.orientation}

  mkdir -p sim_meta batches
  cp sim/simulate.toml sim/loci.tsv sim/regions.tsv sim/tiling.tsv sim_meta/

  # Distribute regions round-robin across N_BATCHES for balanced shard sizes (portable across Linux/macOS)
  tail -n +2 sim/regions.tsv | cut -f1 \
      | awk -v B="${params.n_batches}" '{b = sprintf("%03d", (NR - 1) % B); print > ("batch_ids_" b)}'

  for L in $(echo "${params.read_lengths}" | tr ',' ' '); do
      for bfile in batch_ids_*; do
          bid="${bfile#batch_ids_}"
          awk -v d="sim/fastq/L${L}" '{print d "/" $0 "_R1.fq.gz"}' "$bfile" | xargs cat > "batches/batch_L${L}_${bid}_R1.fq.gz"
          awk -v d="sim/fastq/L${L}" '{print d "/" $0 "_R2.fq.gz"}' "$bfile" | xargs cat > "batches/batch_L${L}_${bid}_R2.fq.gz"
      done
  done

  # Delete the 49,024 individual region FASTQs before Nextflow stages outputs to S3
  rm -rf sim/fastq
  ```

#### Process 3: `ALIGN_BATCH`
- **Resources**: `cpus 16`, `memory '16 GB'`, `disk '20 GB'`, `time '2h'` (Spot-friendly; `errorStrategy 'retry'`, `maxRetries 3`)
- **Inputs**: `tuple val(L), val(batch_id), path(r1), path(r2)` + `path bt2_index_dir`
- **Outputs**: `tuple val(L), path("shard_L${L}_${batch_id}.bam")`
- **Script**:
  ```bash
  bowtie2 -x "${bt2_index_dir}/GRCz12tu" \
      --threads ${task.cpus} \
      -X 500 \
      -1 "$r1" \
      -2 "$r2" \
      2> "shard_L${L}_${batch_id}.bowtie2.log" \
      | samtools view -@ 2 -b -o "shard_L${L}_${batch_id}.bam" -
  ```

#### Process 4: `CAT_AND_SCORE`
- **Resources**: `cpus 4`, `memory '16 GB'`, `disk '160 GB'`, `time '6h'`
  *(Note: `unimap score` threads over read lengths `L=60,150`, so it uses 2 active worker threads and ~8–12 GB RAM total).*
- **Inputs**:
  - `path meta_files` (staged into `sim/`)
  - `path shard_bams` (all shard BAMs for `L=60` and `L=150`)
  - `path fasta`, `path fai`
- **Outputs**:
  - `path "sim/scores/mappability.zarr.tar"`
  - `path "sim/scores/*.parquet"` (`misplacement_summary.parquet`, `bleed_matrix.parquet`, `misplaced_sites.parquet`)
  - `path "sim/scores/score.toml"`
- **Script**:
  ```bash
  mkdir -p sim/aln
  cp ${meta_files} sim/
  bowtie2 --version > sim/aln/bowtie2_version.txt

  for L in $(echo "${params.read_lengths}" | tr ',' ' '); do
      samtools cat -o "sim/aln/sim_L${L}.bam" shard_L${L}_*.bam
      rm -f shard_L${L}_*.bam
  done

  julia -t ${task.cpus} --project=/opt/UniMap.jl -m UniMap score \
      -s sim \
      -f "$fasta" \
      --cutoffs "${params.cutoffs}"

  tar -cf sim/scores/mappability.zarr.tar -C sim/scores mappability.zarr
  rm -rf sim/scores/mappability.zarr
  ```

---

## 5. Container & Environment Requirements

- **Option A — Single Container (simplest)**:
  - Base image with `julia >= 1.12`, `bowtie2` (2.5.x), `samtools`, `awk`, `gzip`, `tar`.
  - Copy this repository (`UniMap.jl`) into `/opt/UniMap.jl` and run:
    ```sh
    julia --project=/opt/UniMap.jl -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
    ```
    at image build time (avoiding package compilation in Batch jobs; note that `GLMakie` is lazy-loaded only by `unimap view` and is not loaded by CLI batch commands).
- **Option B — Per-Process Containers**:
  - `BUILD_LOCI`, `SIMULATE_AND_BATCH`, `CAT_AND_SCORE`: Julia + `UniMap.jl` + `samtools` container.
  - `ALIGN_BATCH`: Standard BioContainers / Wave image with `bowtie2` + `samtools`.

---

## 6. Expected Runtime, Storage & AWS Cost

| Process | Tasks | vCPUs / Task | RAM / Task | Wall Time | Est. AWS Spot Cost |
|---|---|---|---|---|---|
| `BUILD_LOCI` | 1 | 2 | 4 GB | ~30 sec | < $0.01 |
| `SIMULATE_AND_BATCH` | 1 | 32 | 32 GB | ~15–20 min | ~$0.30 |
| `ALIGN_BATCH` | 80 | 16 | 16 GB | ~25–35 min | ~$18–$20 |
| `CAT_AND_SCORE` | 1 | 4 | 16 GB | ~2.5–3 hours | ~$0.40 |
| **Total** | **83** | | | **~3.5–4 hours** | **~$19–$22** |

### Intermediate S3 Storage Footprint
- **Batched FASTQs** (`80 × 2` files): **~65 GB**
- **Shard BAMs** (`80` files): **~110 GB**
- **Final published outputs** (`mappability.zarr.tar` + `.parquet` + `score.toml`): **~7 GB**
