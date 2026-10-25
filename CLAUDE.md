# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Nextflow (DSL2) pipeline that computes empirical paired-end mappability of genomic loci (default: every protein-coding gene in a GFF3) with [UniMap.jl](https://github.com/tdw-89/UniMap.jl). The design point is sharding the bowtie2 step — nearly all of the compute — across many independent tasks while producing results **bit-identical** to UniMap's own monolithic `simulate → align-script → score` path. UniMap.jl is expected as a sibling checkout at `../UniMap.jl` (see `mappability.code-workspace`).

## Commands

Nextflow 25.04 needs Java 17–24. On a Mac with a newer default JDK: `export JAVA_HOME=$(/usr/libexec/java_home -v 24)` (`test/run_tests.sh` does this itself).

```sh
test/run_tests.sh            # planner + check_score.jl unit tests, launch-validation checks, stub run (seconds)
test/run_tests.sh --e2e      # + full fixture run compared exactly against monolithic UniMap (~3 min; needs bowtie2, samtools)
test/run_tests.sh --docker   # + -profile test,docker (slow under amd64 emulation on Apple Silicon)

nextflow run . -profile test,local --unimap_project ../UniMap.jl   # fixture, tools from PATH
nextflow run . -profile test,docker                                # fixture, containers
nextflow run . -profile test,local --unimap_project ../UniMap.jl -stub-run   # wiring only
nextflow run . --help

nextflow config . -profile unity | sed -n '/ALIGN_BATCH/,/}/p'     # inspect effective resources

test/make_fixture.py test/data       # regenerate the fixture deterministically
containers/check_containers.sh       # assert the commands the pipeline runs exist in both images
```

There is no way to run a single test case; `run_tests.sh` is one script of sections. UniMap is found via `$UNIMAP_PROJECT`, default `../UniMap.jl`. Tests run Nextflow from a temp dir so no `.nextflow/` or results land in the repo.

## Architecture

`main.nf` does launch-time validation (`validateParams`, `checkBowtie2Args`, `bowtie2IndexFiles`) and the channel wiring; each process lives in `modules/local/*.nf` with a header comment explaining its invariants and checks — read it before editing a module.

```
BUILD_LOCI → SIMULATE_AND_BATCH ─┬─ fastq batches → ALIGN_BATCH × (n_batches × read lengths) ─┐
                                 └─ sim_meta/ + batches.tsv ─────────────────────────────────→ SCORE
```
(`SAMTOOLS_FAIDX` / `BOWTIE2_BUILD` run only when no `.fai` / `--bowtie2_index` is supplied.)

- **SIMULATE_AND_BATCH** runs `unimap simulate`, then `bin/plan_batches.sh` packs regions into ≤ `n_batches` batches (greedy longest-first by pair count, deterministic). Per-region FASTQs are concatenated into `batch_L{L}_{batch}_R{1,2}.fq.gz` inside the task and never leave it (tens of thousands of small files would dominate an S3 work dir). `batches.tsv` records the pair count every batch must hold.
- `main.nf` parses batch filenames into `[[L, batch], [R1, R2]]` and joins against `batches.tsv` with `failOnMismatch`/`failOnDuplicate`.
- **ALIGN_BATCH** pipes bowtie2 straight into `samtools view -b`, writing `sim_L{L}_{batch}.bam`, and checks bowtie2's input pair count against the plan.
- **SCORE** stages the shards into `sim/aln/` unconcatenated — `unimap score` reads `sim_L{L}_*.bam` as shards and scans them in parallel (needs UniMap with shard support; see `UNIMAP_REF`). It checks every planned shard arrived and all came from one bowtie2 build, writes `align.toml`, runs `unimap score`, then `bin/check_score.jl` verifies every simulated pair was scored exactly once, and finally tars the Zarr store.
- `bin/` is on the task PATH (Nextflow convention): `run_unimap` wraps `julia --project=$UNIMAP_PROJECT -m UniMap`, with threads from `UNIMAP_THREADS`.

### Invariants that must not be broken

A wrong answer here does not look like an error — a pair missing from a BAM is scored as a fragment that failed to map. So:

- **Never sort, filter, or `samtools merge` the BAMs.** `unimap score` streams pair by pair and needs every pair (including unaligned) present with mates adjacent. Shards are scored as they are; there is no need to concatenate them.
- **`--no-unal`, `-X/--maxins`, threads, index, inputs/outputs are forbidden in `--bowtie2_args`** (`checkBowtie2Args`). `-X` is always set to the largest fragment length.
- Sharding is exact only because bowtie2 seeds each pair from its own name/sequence/quality and every read name carries its origin region. Shard order doesn't matter; shard completeness does. Keep every count check (gene IDs → loci, FASTQs → batches, bowtie2 pairs → plan, shards → plan, pairs → score.toml).
- Genes are selected by GFF3 `ID=`, not `Name=` (gene symbols repeat in GRCz12tu). Every selected ID must become a locus or the run stops.
- All processes run under `bash -euo pipefail` (`conf/base.config`); failures upstream of a pipe must not leave a short file.
- Every process has a `stub:` block; keep them in sync when changing outputs so `-stub-run` keeps wiring every process.

### Configuration

- `nextflow.config` holds params (documented inline) and includes `conf/base.config` (resources + retry policy: retry only on OOM/segfault/SIGTERM/Spot-reclaim exit codes, memory doubles per attempt).
- Profiles: `local` (no containers, needs `--unimap_project`), `docker`/`apptainer` (engine only; combine with `local` or `test`), `test` (fixture + scaled grid), `unity` (Slurm + Apptainer), `aws` (direct AWS Batch, needs `--aws_queue`), `seqera`.
- **A `withName` block in a profile replaces the one in `conf/base.config` rather than merging**, so profile blocks restate cpus, memory and time together.
- **`seqera` must carry only resources** — no executor, queue, region, docker, wave, fusion or scratch settings; Seqera's compute environment injects those. Use `aws` only for launches outside Seqera.
- AWS memory requests sit ~4–6% under nominal instance RAM (e.g. 15 GB, not 16) so they fit ECS-available memory.
- Inputs may be `s3://` URIs; always resolve paths with `file()`, never by joining onto the launch dir.

### Containers

- `ghcr.io/tdw-89/unimap:0.2.0` (`containers/unimap.Dockerfile`): Julia + UniMap.jl at a pinned commit (`UNIMAP_REF`), precompiled for a multi-CPU target list, plus samtools. Bump `UNIMAP_REF` and the image tag together; always build `--platform linux/amd64`. Keep the image names in `containers/check_containers.sh` in step with `params.container_*`.
- The bowtie2 biocontainer (bowtie2 2.5.2) differs from Homebrew's bowtie2 used by `-profile local`; results to be compared should come from the same bowtie2 (`score.toml` records which ran).
- Under amd64 emulation on Apple Silicon, multithreaded Julia segfaults — pass `--julia_threads 1`.

### Test fixture

`test/data` (~120 kb, from `test/make_fixture.py`) has planted features — a paralog, an N block, merging loci, a chromosome-edge gene, two genes sharing `Name=`, and a lncRNA and pseudogene that must not be selected. `test/check_scores.jl expect` asserts each one; `check_scores.jl equal` does the sharded-vs-monolithic comparison. If you change the fixture, update both.

## Other files

- `nextflow_aws_proposal.md` — the original design proposal (grid choices, cost estimates).
- `notes.md` — a review checklist; most of it names processes (`BAM_CLEAN`, `BAM_MERGE`, `readBamSheet`, …) that are not in this repo, so check which pipeline an item applies to before acting on it.
