# mappability

Empirical paired-end mappability of a set of loci (by default every protein-coding
gene) with [UniMap.jl](https://github.com/tdw-89/UniMap.jl), as a Nextflow pipeline that
shards the bowtie2 step — nearly all of the compute — across many independent tasks.
It works for any genome: give a FASTA, a GFF3 that names the same sequences, and the read
and fragment lengths of the libraries the mappability is for.

```
GFF + FASTA
     │
 BUILD_LOCI           gene IDs (GFF3 ID=) -> loci.tsv           unimap loci
     │
 SIMULATE_AND_BATCH   tile loci -> per-region FASTQs, packed    unimap simulate
     │                into n_batches balanced batches per L
     ├── simulation metadata ─────────────────┐
     │                                        │
 ALIGN_BATCH × (n_batches × read lengths)     │                 bowtie2 | samtools view
     │   unsorted, unfiltered shard BAMs      │
     │                                        │
 SCORE ◄──────────────────────────────────────┘                 unimap score
         shards, scanned in parallel -> mappability.zarr(.tar) + Parquet summaries
```

Sharding is exact: bowtie2 aligns each pair independently (its tie-breaking PRNG is
seeded from the read's own name, sequence and quality), and every read name carries the
region it was simulated from, so the shards, in any order, score identically to one
monolithic run. `unimap score` reads the shards directly (`sim_L{L}_{batch}.bam`) and scans
them in parallel, so they are never concatenated. `test/run_tests.sh --e2e` checks this: every Zarr array and Parquet
table from the pipeline is compared with UniMap's own `simulate → align-script → score` path.
What must never happen is sorting or filtering the BAMs: `unimap score` streams them pair by
pair and needs every pair present with its mates adjacent.

## Quick start

```sh
# The bundled mini-genome, no containers (julia, bowtie2, samtools on PATH):
nextflow run . -profile test,local --unimap_project ../UniMap.jl

# The same in containers:
nextflow run . -profile test,docker

# Real data on a workstation -- a handful of genes, not the genome:
nextflow run . -profile local --unimap_project ../UniMap.jl \
    --fasta ref/genome.fa --gff ref/genes.gff3 --bowtie2_index ref/bowtie2/genome \
    --read_lengths 50,150 --fragment_lengths 100,150,200,300,400,500 \
    --genes my_genes.txt --n_batches 4

# Every protein-coding gene, on Unity:
nextflow run . -profile unity \
    --fasta ... --gff ... --bowtie2_index ... \
    --read_lengths ... --fragment_lengths ... \
    -w /scratch/workspace/$USER-mappability/work --outdir results
```

The zebrafish run this was built for (GRCz12tu, NCBI RefSeq annotation, 28,407 genes) used
`--read_lengths 60,150 --fragment_lengths 50,75,100,125,150,175,200,250,300,350,400,450,500`.

`nextflow run . --help` lists every parameter. Nextflow 25.04 needs Java 17–24; on a Mac
with a newer default JDK, `export JAVA_HOME=$(/usr/libexec/java_home -v 24)`.

## Parameters

| Parameter | Default | |
|---|---|---|
| `--fasta` | required | uncompressed reference; a `.fai` beside it is reused, otherwise built |
| `--gff` | required | GFF3 from any source (or GTF with `--genes`), optionally `.gz` |
| `--read_lengths` | required | comma-separated, e.g. `50,150` |
| `--fragment_lengths` | required | comma-separated; every (L, F) is simulated, F < L allowed; bowtie2 `-X` is set to the largest |
| `--bowtie2_index` | build one | index prefix built from `--fasta`; omit only for small genomes |
| `--gene_biotype` | `protein_coding` | select every GFF3 gene of this biotype (comma-separate several), by `ID=` |
| `--biotype_attribute` | `auto` | `gene_biotype` (NCBI), `biotype` (Ensembl) or `gene_type` (GENCODE), whichever genes carry; with none, protein-coding = has a CDS |
| `--genes` | — | or a file of IDs, one per line (anything `unimap loci -g` accepts) |
| `--upstream` / `--downstream` | `1500` / `1500` | strand-aware flanks |
| `--step` | `10` | bp between fragment starts |
| `--orientation` | `plus` | `plus`, `minus` or `both` |
| `--pairs_per_batch` | `35000000` | read pairs per alignment task (~1 h of bowtie2 on 16 cores); sets the number of batches |
| `--n_batches` | — | or a fixed number of batches per read length (≤ 999, capped at the region count) |
| `--bowtie2_args` | — | extra bowtie2 arguments; `--no-unal`, `-X`, threads, I/O are refused |
| `--cutoffs` | `0,2,5,10,20,30` | fragment MAPQ cutoffs |
| `--frag_mapq` / `--tolerance` | `min` / `0` | passed to `unimap score` |
| `--tar_zarr` | `true` | publish `mappability.zarr` as one tar (no recompression) |
| `--publish_bams` | `false` | also publish the shard BAMs (~110 GB at genome scale) |
| `--outdir` | `results` | local path or `s3://` URI |
| `--unimap_project` | `$UNIMAP_PROJECT` | UniMap.jl checkout; only without containers |
| `--julia_threads` | task cpus | set `1` only to run the amd64 image under emulation (Apple Silicon), where multithreaded Julia segfaults |

### Which genes

GFF3 has no standard biotype attribute, so `bin/select_genes.sh` looks at what the genes
(root features with an `ID=` whose type contains `gene`) carry:

| Annotation | Attribute | Gene IDs |
|---|---|---|
| NCBI RefSeq | `gene_biotype=protein_coding` | `gene-…` |
| Ensembl | `biotype=protein_coding` | `gene:ENSG…` |
| GENCODE | `gene_type=protein_coding` | `ENSG….N` |
| none (MAKER, BRAKER, FlyBase-style, …) | — | a gene with a CDS (not `pseudo=true`, not a `pseudogene`) is protein-coding |

`loci/gene_selection.tsv` records the rule that was used. Genes on sequences the FASTA
doesn't have (alt haplotypes, patches, scaffolds left out of a primary-assembly FASTA) can't be
simulated. They are listed in `loci/excluded_genes.tsv` with a warning. If that leaves no genes, the
GFF and FASTA name their sequences differently (`chr1` vs `1` vs `NC_…`) and the run stops.
Sequence names may not contain `|` (UniMap's read names are `|`-separated).

Genes are selected by `ID=`, not `Name=`: symbols repeat (23 protein-coding symbols occur
twice in zebrafish GRCz12tu), and `unimap loci` skips ambiguous aliases. Every selected ID
must become a locus or the run stops.

## Output

```
results/
  loci/          gene_ids.txt, gene_selection.tsv, excluded_genes.tsv, loci.tsv
  simulation/    simulate.toml, loci.tsv, loci.bed, regions.tsv, tiling.tsv,
                 batch_plan.tsv (region -> batch), batches.tsv (pairs per batch)
  alignment/logs/sim_L{L}_{batch}.bowtie2.log
  alignment/bams/sim_L{L}_{batch}.bam                (--publish_bams only)
  scores/        mappability.zarr.tar, misplacement_summary.parquet,
                 bleed_matrix.parquet, misplaced_sites.parquet,
                 score.toml, align.toml
  pipeline_info/ timeline, report, trace, DAG
```

`tar -xf mappability.zarr.tar` gives the store `unimap view` and `locus_track` read.
`score.toml` records the UniMap commit, reference checksum, bowtie2 version and — through
`align.toml` — the aligner arguments and the Nextflow run that produced it.

## Checks

A wrong answer here does not look like an error: a pair missing from a BAM is scored as a
fragment that failed to map. So the pipeline counts:

- **BUILD_LOCI** — every selected gene ID became a locus.
- **SIMULATE_AND_BATCH** — every per-region FASTQ was packed into exactly one batch.
- **ALIGN_BATCH** — bowtie2's input pair count equals the planner's count for that batch.
- **SCORE** — every planned shard arrived; all were aligned by one bowtie2 build;
  and after scoring, every simulated pair was seen exactly once with both mates
  (`bin/check_score.jl` against `score.toml`).

## Profiles

| Profile | |
|---|---|
| `local` | workstation, no containers; tools from `PATH`, UniMap from `--unimap_project` |
| `docker` / `apptainer` | container engine; combine with `local` or `test` |
| `test` | bundled ~120 kb fixture (`test/data`), scaled-down grid |
| `unity` | Slurm + Apptainer on UMass Unity |
| `aws` | AWS Batch launched directly with Nextflow; needs `--aws_queue` |
| `seqera` | Seqera Platform: resources only (see below) |

### Seqera Platform (AWS Batch)

Select the **`seqera`** profile, not `aws`. Seqera's compute environment supplies the
executor, queue, region, work directory, Wave and Fusion; `seqera` carries only cpus, memory
and time, plus local scratch for two tasks (`conf/fusion_scratch.config`). Recommended
compute environment: Fusion + Wave + fast instance storage, Spot, instance families with
local NVMe (`c6id`, `m6id`, `r6id`).

Under Fusion a task runs inside its S3 work directory, so every file it writes is an S3
object. `SIMULATE_AND_BATCH` (49k per-region FASTQs) and `SCORE` (the ~450k-file Zarr store
before it is tarred) therefore run in local scratch — the NVMe, mounted as the container's
`/tmp` — and only their declared outputs go to S3: ~70 GB of local disk for
`SIMULATE_AND_BATCH`, ~15 GB for `SCORE` (the shards stream from S3). Both log `df -h .` to
`.command.err`. Without fast instance storage, scratch is a 100 GB boot disk, too small at
genome scale. Without Fusion every task downloads its inputs first, and `SCORE` needs
~125 GB of EBS — see `conf/aws.config`.

Every input may be an `s3://` URI. `--outdir` on S3 publishes by server-side copy.

### Unity

Export `APPTAINER_TMPDIR`, `APPTAINER_CACHEDIR` and `NXF_APPTAINER_CACHEDIR` to `/work`
directories before launching, and put `-w` on scratch: at genome scale the work directory
holds ~65 GB of batched FASTQ and ~110 GB of shard BAMs.

## Containers

| Image | Used by | |
|---|---|---|
| `ghcr.io/tdw-89/unimap:0.2.0` | loci, simulate, score, faidx | `containers/unimap.Dockerfile`: Julia 1.13.2 + UniMap.jl (pinned commit) precompiled, samtools |
| nf-core's `mulled-v2-ac74a7…:f70b31…-0` | bowtie2 build/align | bowtie2 2.5.2, samtools 1.18 |

Build and push the UniMap image (always `linux/amd64` — Unity and Batch are x86_64):

```sh
docker buildx build --platform linux/amd64 -f containers/unimap.Dockerfile \
    -t ghcr.io/tdw-89/unimap:0.2.0 --push .
containers/check_containers.sh          # or: sbatch containers/check_containers.sh
```

Precompiled caches are built for the multi-target CPU list of the official Julia binaries,
so they load on any Batch instance or Unity node instead of recompiling in every task.

bowtie2 versions can place reads differently. The `local` profile uses whatever bowtie2 is
on `PATH` (2.5.5 from Homebrew), the containers use 2.5.2; `score.toml` records which one
ran, and results to be compared should come from the same one.

## Tests

```sh
test/run_tests.sh            # planner, score check, launch validation, stub run
test/run_tests.sh --e2e      # + full fixture run, compared exactly against monolithic UniMap
test/run_tests.sh --docker   # + the containerised run
```

`test/make_fixture.py` regenerates `test/data` deterministically. The fixture has a planted
paralog, an N block, loci whose flanks merge into one region, a gene at a chromosome edge,
two genes sharing a `Name=`, and a lncRNA and a pseudogene that must not be selected;
`test/check_scores.jl expect` checks each comes out as it must.

## Scale

The work is set by the read pairs, roughly per read length

    pairs ≈ (bp of all regions) × (number of fragment lengths) × (2 if orientation=both) / step

and the run log prints the actual count once `SIMULATE_AND_BATCH` is done. Alignment is
sharded by `--pairs_per_batch`, so each `ALIGN_BATCH` task takes about the same time whatever
the genome. Its memory is mostly the bowtie2 index (~2.5 GB per Gb of genome), so 28 GB
covers genomes up to ~8 Gb. The serial steps do grow with the pair count:
`SIMULATE_AND_BATCH` disk (~25 bytes per pair), and `SCORE` memory and time. The default
resources were sized for the zebrafish run's ~2.8 billion pairs in all. For much more than that,
raise them with a `-c` config (restate cpus, memory and time in each `withName` block).

### The zebrafish run, for reference

28,407 genes → 12,256 regions (1,060 Mb); ~1.38 billion pairs per read length; ~1,000 vCPU-hours of
bowtie2. The wall time is two serial steps around a parallel one: `SIMULATE_AND_BATCH`
(~20 min on 32 cores), 80 × `ALIGN_BATCH` (40 batches × 2 read lengths), then `SCORE`. Measured components of `SCORE`,
extrapolated to genome scale on an Apple M-series laptop: the shard scans and per-region
tracks take a few minutes on 8 threads; writing the Zarr store takes ~12 min per read length
and is serial (Zarr.jl's Blosc codec uses global state), so ~25–30 min in all. The store holds
~600 GB of counts per read length before compression, and copying that into chunks, not the
compression, is most of the cost.
