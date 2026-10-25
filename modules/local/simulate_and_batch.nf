/*
 * Tile every locus with error-free read pairs, then pack the per-region FASTQs
 * into a few large batches.
 *
 * `unimap simulate` writes one R1/R2 pair per simulation region per read length
 * -- 49,024 files for the zebrafish grid. Staged through an S3 work directory
 * that many small objects would dominate the run, so they never leave this
 * task: they are concatenated into batch_L{L}_{batch}_R{1,2}.fq.gz (gzip
 * members concatenate into a valid gzip stream) and deleted as they are
 * consumed, which also keeps peak disk at roughly one copy of the reads.
 * Under Fusion even the work directory is S3, so conf/fusion_scratch.config runs
 * this task in local scratch; the `df` line in .command.err shows where.
 *
 * R1 and R2 of a batch are built from the same region list in the same order,
 * so mates stay in step. Which batch a region lands in, and in what order, does
 * not matter to scoring: every read name carries its region, and bowtie2 seeds
 * each pair from its own name and sequence, so sharded alignment is identical
 * to one monolithic run.
 *
 * sim_meta/ is all `unimap score` needs from here -- it re-enumerates the
 * fragments from the reference rather than reading the FASTQs.
 */

process SIMULATE_AND_BATCH {
    tag   { "L=${params.read_lengths}" }

    container params.container_unimap

    publishDir "${params.outdir}/simulation", mode: params.publish_mode,
               saveAs: { fn -> fn.startsWith('sim_meta/') ? fn - 'sim_meta/' : null }

    input:
    path fasta
    path fai
    path loci

    output:
    path "sim_meta/*",                emit: meta
    path "sim_meta/batches.tsv",      emit: manifest
    path "batches/batch_L*.fq.gz",    emit: fastq

    script:
    def project = params.unimap_project ? "export UNIMAP_PROJECT='${file(params.unimap_project)}'" : ''
    def Ls = params.read_lengths.toString().tokenize(',').collect { it.trim() }.join(' ')
    """
    ${project}
    export UNIMAP_THREADS=${params.julia_threads ?: task.cpus}

    # Where the ~70 GB of FASTQ will land (local scratch, not S3, under Fusion).
    echo "task dir: \${PWD}" >&2
    df -h . >&2 || true

    run_unimap simulate \\
        -f "${fasta}" \\
        -l "${loci}" \\
        -o sim \\
        -L "${params.read_lengths}" \\
        -F "${params.fragment_lengths}" \\
        --step ${params.step} \\
        --orientation ${params.orientation}

    mkdir -p sim_meta batches
    cp sim/simulate.toml sim/loci.tsv sim/loci.bed sim/regions.tsv sim/tiling.tsv sim_meta/

    plan_batches.sh sim/regions.tsv sim/tiling.tsv ${params.n_batches ?: 0} ${params.pairs_per_batch} \\
        sim_meta/batch_plan.tsv sim_meta/batches.tsv

    for L in ${Ls}; do
        for batch in \$(tail -n +2 sim_meta/batch_plan.tsv | cut -f2 | sort -u); do
            awk -F'\\t' -v b="\${batch}" 'NR > 1 && \$2 == b { print \$1 }' sim_meta/batch_plan.tsv > ids.txt
            for mate in R1 R2; do
                # NUL-separated, so no character in a region ID is special to xargs.
                awk -v d="sim/fastq/L\${L}/" -v m="_\${mate}.fq.gz" '{ print d \$0 m }' ids.txt | tr '\\n' '\\0' > files.txt
                xargs -0 cat < files.txt > "batches/batch_L\${L}_\${batch}_\${mate}.fq.gz"
                xargs -0 rm -f < files.txt
            done
        done
    done
    rm -f ids.txt files.txt

    # Every region file must have been consumed by exactly one batch.
    leftover=\$(find sim/fastq -type f | wc -l | tr -d ' ')
    if [ "\${leftover}" -ne 0 ]; then
        echo "error: \${leftover} FASTQ file(s) were not assigned to any batch:" >&2
        find sim/fastq -type f | head >&2
        exit 1
    fi
    rm -rf sim/fastq
    """

    stub:
    def Ls = params.read_lengths.toString().tokenize(',').collect { it.trim() }.join(' ')
    """
    mkdir -p sim_meta batches
    touch sim_meta/simulate.toml sim_meta/loci.tsv sim_meta/regions.tsv sim_meta/tiling.tsv
    printf 'L\\tbatch\\tn_regions\\tn_pairs\\n' > sim_meta/batches.tsv
    for L in ${Ls}; do
        printf '%s\\t001\\t1\\t0\\n' "\${L}" >> sim_meta/batches.tsv
        touch batches/batch_L\${L}_001_R1.fq.gz batches/batch_L\${L}_001_R2.fq.gz
    done
    """
}
