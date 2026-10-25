## Please check in the code

### 1. `-profile aws` versus the Seqera compute environment

Read `nextflow.config` and `conf/aws.config`. Seqera already injects the following from the compute environment:

- executor (`awsbatch`)
- queue
- region
- work dir
- Wave/Fusion settings

Report what the `aws` profile sets and whether any of it conflicts with that injected config. In particular, check these settings:

- `process.queue` / `aws_cpu_queue`
- `aws.region`
- `aws.batch.cliPath` (not applicable to Batch Forge + Fusion instances)
- `docker.*`, `wave.*`, `fusion.*`
- `process.scratch`

Then recommend one of:

- **(a)** select no profile in Seqera, or
- **(b)** add a minimal `seqera` profile that carries only pipeline-specific settings (resources, retries, labels) and leaves executor, queue and storage to the compute environment.

### 2. `readBamSheet` with `s3://` URIs

The README says relative paths resolve against the launch directory and every path is checked at launch. Confirm the following:

- An `s3://bucket/key` string is treated as absolute and never joined onto `launchDir`. An absoluteness test like `startsWith('/')` would break it.
- The existence check works on S3 through Nextflow's `file()`.
- The `--bam_pattern` warning matches on the object name.
- `--merge_input` itself can be an `s3://` path.

Add a test case if it's cheap.

### 3. Symlinks and link-mode publishing on an object store

`BAM_CLEAN` passes an input with nothing to remove through as a **symlink**, to avoid rewriting tens of GB. On an S3 work dir, check what actually happens, both with and without Fusion. The possibilities are:

- the link survives,
- it is silently dereferenced into a full copy (tens of GB per input), or
- it breaks `BAM_MERGE`.

Also check `--merged_bam_publish_mode link` with an S3 `--outdir`. It should probably be forced to, or documented as, `copy` there. On S3, `copy` is a server-side copy.

### 4. Disk, memory and temp space on Batch

**Without Fusion,** each task downloads all its inputs to local disk before running and uploads outputs afterwards. That means:

- `BAM_MERGE` would need every cleaned input plus the merged output on local disk at once.
- Every QC task (cramino, mosdepth, bam_scan, modkit summary) and `MODKIT_PILEUP` would each download the whole merged BAM.

Estimate the peak local disk per process for this dataset, and say whether the non-Fusion compute environment could run it at all.

**`BAM_READ_NAMES`:** BusyBox sort runs in memory, at about 2.5 GB per 10M reads (sized in `conf/base.config`). If the request is computed from input size at launch, confirm that works for S3 paths. Also confirm the largest request fits the compute environment's instance types.

**`CROSS_INPUT_OVERLAP`:** GNU sort spills to disk. Find where it spills (`-T` / `TMPDIR`). Spilling into the Fusion-mounted S3 work dir would be slow; local NVMe or `/tmp` is what's wanted.

### 5. Containers under Wave

Fusion requires Wave, which adds the Fusion layer to each task image. Confirm the pinned biocontainers (BusyBox samtools and the python image) run under it, and that no script depends on a host path that exists on Unity but not on Batch.
