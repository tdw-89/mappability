/*
 * Score every simulated fragment from the alignment shards.
 *
 * The shards go to `unimap score` as they are -- sim/aln/sim_L{L}_{batch}.bam --
 * and it scans them in parallel. Nothing is concatenated: each shard holds whole
 * pairs with their mates adjacent, exactly as bowtie2 wrote them, and shard
 * order is irrelevant because every read name carries the region it came from.
 * (UniMap refuses a pair that turns up in two shards.)
 *
 * `-f` is passed explicitly: simulate.toml records the FASTA's absolute path
 * inside the SIMULATE_AND_BATCH task, which does not exist here.
 *
 * Checks, in order:
 *   - every batch in the plan has its shard (a missing one would otherwise
 *     score as hundreds of thousands of unmappable fragments)
 *   - every shard was aligned by the same bowtie2 build
 *   - after scoring, every simulated pair was seen exactly once
 *     (bin/check_score.jl against score.toml)
 *
 * align.toml is written into the simulation directory so score.toml records
 * the aligner settings, as it would after `unimap align-script`.
 *
 * Disk: only what this task writes -- the Zarr store and its tar, ~7 GB each at
 * genome scale. Under Fusion the shards are read straight from S3; without it
 * (plain AWS Batch) they are downloaded first, ~110 GB. The unpacked store is
 * ~450k chunk files, so under Fusion conf/fusion_scratch.config runs this task
 * in local scratch and only the tar goes to S3; the `df` line in .command.err
 * shows where.
 */

process SCORE {
    tag   { "L=${params.read_lengths}" }

    container params.container_unimap

    publishDir "${params.outdir}/scores", mode: params.publish_mode,
               saveAs: { fn -> fn.tokenize('/').last() }

    input:
    path meta,     stageAs: 'sim/*'
    path manifest, stageAs: 'batches.tsv'
    path shards,   stageAs: 'sim/aln/*'
    path versions, stageAs: 'versions/bowtie2_version*.txt'
    path fasta
    path fai
    val  index_name
    val  aligner_args

    output:
    path "sim/scores/mappability.zarr${params.tar_zarr ? '.tar' : ''}", emit: zarr
    path "sim/scores/*.parquet",                                       emit: tables
    path "sim/scores/score.toml",                                      emit: toml
    path "sim/align.toml",                                             emit: align_toml

    script:
    def project = params.unimap_project ? "export UNIMAP_PROJECT='${file(params.unimap_project)}'" : ''
    def toml_args = aligner_args.tokenize(' ').collect { '"' + it.replace('\\', '\\\\').replace('"', '\\"') + '"' }.join(', ')
    """
    ${project}
    export UNIMAP_THREADS=${params.julia_threads ?: task.cpus}

    # Where the Zarr store will land (local scratch, not S3, under Fusion).
    echo "task dir: \${PWD}" >&2
    df -h . >&2 || true

    # ---- every planned shard is here --------------------------------------
    missing=0
    while IFS=\$'\\t' read -r L batch n_regions n_pairs; do
        if [ ! -e "sim/aln/sim_L\${L}_\${batch}.bam" ]; then
            echo "error: no shard for L=\${L} batch \${batch} (\${n_pairs} pairs)" >&2
            missing=\$(( missing + 1 ))
        fi
    done < <(tail -n +2 batches.tsv)
    n_planned=\$(( \$(wc -l < batches.tsv) - 1 ))
    n_shards=\$(ls sim/aln/ | grep -c '\\.bam\$' || true)
    if [ "\${missing}" -ne 0 ] || [ "\${n_shards}" -ne "\${n_planned}" ]; then
        echo "error: \${n_planned} batch(es) planned, \${n_shards} shard(s) received" >&2
        exit 1
    fi

    # ---- one bowtie2 build for all of them --------------------------------
    first=\$(ls versions/* | head -n 1)
    for v in versions/*; do
        if ! cmp -s "\${first}" "\${v}"; then
            echo "error: shards were aligned by different bowtie2 builds:" >&2
            cat versions/* | grep -i version | sort | uniq -c >&2
            exit 1
        fi
    done
    cp "\${first}" sim/aln/bowtie2_version.txt

    {
        printf '%s\\n' \\
            'aligner = "bowtie2"' \\
            'index = "${index_name}"' \\
            'aligner_args = [${toml_args.replace("'", "'\\''")}]' \\
            'bowtie2 = "bowtie2"' \\
            'samtools = "samtools"' \\
            'nextflow_session = "${workflow.sessionId}"' \\
            'nextflow_run = "${workflow.runName}"' \\
            'pipeline = "${workflow.manifest.name} ${workflow.manifest.version}"' \\
            'pipeline_commit = "${workflow.commitId ?: ''}"'
        printf 'script = "ALIGN_BATCH: %s shards, scored as shards"\\n' "\${n_planned}"
        printf 'created = "%s"\\n' "\$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } > sim/align.toml

    # ---- score ------------------------------------------------------------
    run_unimap score \\
        -s sim \\
        -f "${fasta}" \\
        --cutoffs "${params.cutoffs}" \\
        --frag-mapq ${params.frag_mapq} \\
        --tolerance ${params.tolerance}

    check_score.jl sim/simulate.toml sim/scores/score.toml

    if ${params.tar_zarr}; then
        tar -cf sim/scores/mappability.zarr.tar -C sim/scores mappability.zarr
        rm -rf sim/scores/mappability.zarr
    fi
    """

    stub:
    """
    mkdir -p sim/scores
    touch sim/scores/mappability.zarr${params.tar_zarr ? '.tar' : ''} sim/scores/misplacement_summary.parquet \\
          sim/scores/score.toml sim/align.toml
    """
}
