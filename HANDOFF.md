# Handoff: mappability pipeline

Status as of 2026-10-24, for whoever picks this up next. `README.md` (usage) and `CLAUDE.md`
(architecture, invariants) are current; this file covers state, findings and next steps.

## What this is

A Nextflow pipeline (`main.nf`, `modules/local/`, `conf/`) that runs
[UniMap.jl](https://github.com/tdw-89/UniMap.jl) (sibling checkout `../UniMap.jl`) to compute
empirical paired-end mappability for all 28,407 zebrafish GRCz12tu protein-coding genes, as
specified in `nextflow_aws_proposal.md`. The bowtie2 step (~1,000 vCPU-h) is sharded into
`n_batches` × read-length tasks. Targets: AWS Batch via Seqera Platform (`-profile seqera`)
and UMass Unity (`-profile unity`); `-profile local` runs on a Mac without containers.

```
BUILD_LOCI → SIMULATE_AND_BATCH ─┬─ batch FASTQs → ALIGN_BATCH × (batches × L) ─┐
                                 └─ sim_meta/ + batches.tsv ────────────────────→ SCORE
```

## State

| Item | State |
|---|---|
| UniMap.jl | `3959a74`, **pushed**: `unimap score` reads BAM shards `sim_L{L}_*.bam` in parallel; Blosc-threaded Zarr writes |
| Image `ghcr.io/tdw-89/unimap:0.2.0` | Built locally (amd64), pins `3959a74`, passes `containers/check_containers.sh`. **Not pushed yet**: the user pushes it (`docker push ghcr.io/tdw-89/unimap:0.2.0`) |
| Image `ghcr.io/tdw-89/unimap:0.1.0` | Pushed; the old UniMap (`8db600a`), can't read shards. The pipeline no longer works with it |
| mappability repo | Last commit `52d0c3f` (the user's). The `CAT_AND_SCORE` → `SCORE` refactor, the image bump to 0.2.0, `CLAUDE.md` and this file are **uncommitted** |
| Tests | `test/run_tests.sh --docker`: 64 passed. UniMap `julia -t 4 --project test/runtests.jl`: 238 passed (also at `-t 1`) |
| Real data | Smoke run with 11 random genes, full production grid, `-profile local`: all 813,285 pairs per L scored, 5.5 min. Not yet re-run since the `SCORE` refactor (fixture e2e and a direct UniMap shard-vs-single comparison on those BAMs both match exactly) |
| Seqera / Unity | Never run. The user has a Seqera CE (Fusion on, m5 family, Max CPUs 64) that is **unsuitable**; see below |

Reference data (read-only, the user's): `/Users/thomaswolfe/temp_for_analysis/fish/zebrafish/mappability_broad_marks/ref_links/`
(FASTA + `.fai`, GFF, `bowtie2/GRCz12tu.*.bt2`).

## Verified facts (measured, not assumed)

- The proposal's gene count holds: 28,407 protein-coding IDs → 28,407 loci; `unimap loci` takes 14 s and 1.3 GB RSS.
- **Sharding is exact.** The pipeline's sharded output equals UniMap's monolithic
  `simulate → align-script → align.sh → score` on every Zarr array and Parquet table
  (`test/run_tests.sh --e2e`, `test/check_scores.jl equal`).
- **bowtie2 throughput** is ~600–800 pairs/s per M-series core, which extrapolates to ~1,100 core-h for the full grid (the proposal says ~1,000 vCPU-h).
- **Score cost at genome scale** (extrapolated from profiling):
  - Reading the BAMs is cheap: ~5 min per L single-threaded, now parallel over shards.
  - Writing the Zarr store dominates: ~12 min per L with Blosc threads (21 min before), and it runs serially. Total ~30 min.
  - About 84% of the Zarr time is Zarr.jl copying data into chunk buffers (`copyto!` in `ZarrCore` `writeblock!`); compression is the smaller part. The store holds ~600 GB of uncompressed counts per L (169 counts/position: 13 F + 13×6 UInt16 + 13×6 UInt32).
  - The only big lever left is the storage format (e.g. UInt16 `n_ectopic`, fewer cutoffs). That is a design decision for the user; not done.
- **The Zarr write lock in UniMap is required.** Zarr.jl's Blosc codec calls `blosc_compress` and `blosc_set_compressor`, which use Blosc's global, non-thread-safe state.
- **Under Fusion, task code doesn't run on local disk.** Nextflow runs it inside the Fusion-mounted S3 work dir; `scratch` defaults to false (checked in Nextflow 25.04.6 `FusionScriptLauncher`). So anything written in the work dir goes to S3.
- **Fusion without "fast instance storage"** uses a 100 GB gp3 boot disk whose `/tmp` cache is shared by all tasks on the instance (Seqera docs).

## Next steps, in order

1. **Push the image and commit.** The user pushes `ghcr.io/tdw-89/unimap:0.2.0` and commits the mappability changes (ask first; they commit themselves).
2. **Fix the Fusion scratch issue.** Not done yet; the user deferred it. Two steps write throwaway files into the S3-backed work dir:
   - `SIMULATE_AND_BATCH`: `unimap simulate` writes 49,024 per-region FASTQs under `sim/fastq` before packing them.
   - `SCORE`: `unimap score` writes the unpacked `mappability.zarr` (~450k chunk files) before the tar.

   Proposed fix: set `scratch = true` for those two processes in `conf/seqera.config` (an explicit `scratch` is honored under Fusion). Only declared outputs then get copied to the work dir. Add a `df -h .` line to both scripts to confirm where scratch lands (container `/tmp` may sit on the boot EBS, not NVMe). If it lands on a small disk, use `scratch = '/some/nvme/path'` with a matching `aws.batch.volumes` mount.

   Note `conf/seqera.config`'s comment says that profile must carry "no scratch settings". That rule was written to keep environment settings out; per-process scratch for these two is a pipeline need, so update the comment along with it.

   A deeper alternative: add a UniMap `simulate --n-batches` option that writes batched FASTQs directly. That removes the 49k files and the 65 GB re-read/re-write entirely, but `align-script` would need to understand batches.
3. **Create a new Seqera compute environment.** Don't reuse the user's m5 one: m5 has no NVMe, and Max CPUs 64 runs only 4 alignment tasks at a time (~10–13 h).
   - Batch Forge, in the S3 bucket's region.
   - Wave + Fusion v2 + **fast instance storage**.
   - Fusion snapshots off.
   - Spot, with the head job on demand.
   - Instance types `c6id`, `m6id` (optionally `r6id`).
   - Max CPUs 640–1280; check the account's EC2 Spot vCPU quota first.
   - Allocation strategy `SPOT_PRICE_CAPACITY_OPTIMIZED`.

   Request fit:

   | Process | Request | Fits on |
   |---|---|---|
   | ALIGN_BATCH | 16 CPU / 28 GB | `c6id.4xlarge` |
   | SIMULATE_AND_BATCH | 32 / 60 | `c6id.8xlarge` |
   | SCORE | 8 / 28 | `m6id.2xlarge` |
4. **Roll out in stages:**
   1. `-profile test,seqera` on the new CE.
   2. ~50 genes with `--genes`. Inspect `pipeline_info/trace.txt` (peak RSS, especially SCORE and SIMULATE) and the S3 work-dir object count.
   3. The full run. Inputs (FASTA, `.fai`, GFF, bowtie2 index) are uploaded to S3 in the same region.
5. **Unity:** run `sbatch containers/check_containers.sh` (with the APPTAINER_* variables exported) after pushing 0.2.0, then the test profile.
6. **Re-tune resources from real traces.** The SCORE memory (28–32 GB) is an estimate: fragment tables ~6 GB plus parallel per-region track buffers, which can reach ~1 GB each for >1 Mb regions.

## Gotchas already hit (don't rediscover them)

- **Nextflow 25.04 needs Java 17–24**; the Mac defaults to 26: `export JAVA_HOME=$(/usr/libexec/java_home -v 24)`.
- **Image base must be Debian trixie**, not bookworm: Julia's bundled zstd needs libzstd ≥ 1.5.7 to unpack the registry.
- **Precompile with `Base.compilecache(Base.identify_package("UniMap"))`, not `using UniMap`.**
  - `using UniMap` precompiles all of UniMap's `[deps]`, including GLMakie/CairoMakie. One of those fails headless, and Julia then silently loads UniMap uncached.
  - The Dockerfile's last `RUN` fails the build if the caches aren't reused.
  - Caches are built for the multi-CPU target list `JULIA_CPU_TARGET`.
- **amd64 builds and containers on Apple Silicon run under emulation**, where multithreaded Julia segfaults.
  - The Dockerfile uses one thread and serial downloads at build time.
  - Container tests use `--julia_threads 1`; `test/run_tests.sh` does this automatically on arm64.
  - Not an issue on real x86.
- **bowtie2 versions differ:** the local profile uses Homebrew 2.5.5, the containers 2.5.2 (nf-core's mulled image). Compare results only within one version.
- **`groupTuple(sort: ...)` misbehaved** in this Nextflow; the R1/R2 order is sorted in a following `map` and asserted.
- **Nextflow strips the script's common indentation**, so no heredocs with indented terminators, and no multi-line strings interpolated into `script:` blocks.
- **Julia benchmark scripts need `const` globals.** An untyped global in an inner loop made a Zarr benchmark look ~50× slower than it is.
- **The disk is nearly full** (~35 GB free). Docker build cache from earlier image builds can be pruned (`docker builder prune`), which also clears the user's own cache. Ask first.

## Where things are

- **Pipeline:** `main.nf` (validation + wiring), `modules/local/*.nf` (each has a header explaining its checks), `bin/` (`run_unimap`, `plan_batches.sh` (greedy largest-first batch packing), `check_score.jl`), `conf/*.config`.
- **Tests:** `test/run_tests.sh`, `test/check_scores.jl`, `test/make_fixture.py` → `test/data/` (FASTA and GFF only; the bowtie2 index is built by `BOWTIE2_BUILD`).
- **Containers:** `containers/unimap.Dockerfile`, `containers/check_containers.sh`.
- **UniMap changes:** `src/score.jl` (`bam_paths`, `merge_scans`, the three-phase `score_simulation`), plus a shard testset in `test/runtests.jl`.
- **Benchmark scripts and the smoke-run work dir** are in the previous session's folder (`~/.copilot/session-state/6bb13513-…/files/`). They aren't needed, and the paths there are disposable.
