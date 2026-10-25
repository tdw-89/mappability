/*
 * Align one batch of simulated pairs genome-wide.
 *
 * The BAM is left exactly as bowtie2 writes it: unsorted, unfiltered, every
 * pair present including unaligned ones (no --no-unal), mates adjacent.
 * `unimap score` streams it pair by pair and counts a missing pair as a
 * fragment that failed to map, so a dropped record is a wrong answer, not an
 * error. Sorting would separate mates and make the BAM unscorable.
 *
 * THE CHECK: bowtie2 reports how many pairs it read, and the planner recorded
 * how many this batch's FASTQ holds. They must agree, or the shard is short --
 * a truncated upload, a partial FASTQ -- and every fragment in the gap would
 * score as unmappable.
 *
 * Shards are independent: bowtie2 places each pair on its own, with the
 * tie-breaking PRNG seeded from the read's name, sequence and quality, so
 * splitting the reads changes nothing about any alignment. The sharded-vs-
 * monolithic equivalence test (test/run_tests.sh --e2e) checks exactly that.
 *
 * The shard is named sim_L{L}_{batch}.bam, the form `unimap score` reads as
 * shards of read length L without concatenating them.
 */

process ALIGN_BATCH {
    tag   { "L=${L} batch ${batch}" }

    container params.container_bowtie2

    publishDir "${params.outdir}/alignment/logs", mode: params.publish_mode,
               pattern: "*.bowtie2.log"
    // ~110 GB across all shards at genome scale; only needed to re-score.
    publishDir "${params.outdir}/alignment/bams", mode: params.publish_mode,
               pattern: "*.bam", enabled: params.publish_bams as boolean

    input:
    tuple val(L), val(batch), path(r1), path(r2), val(n_pairs)
    path index_files
    val  index_name
    val  aligner_args

    output:
    path "sim_L${L}_${batch}.bam",          emit: bam
    path "sim_L${L}_${batch}.bowtie2.log",  emit: log
    path "bowtie2_version.txt",               emit: version

    script:
    def prefix = "sim_L${L}_${batch}"
    """
    bowtie2 --version > bowtie2_version.txt

    # bowtie2 gets every core; samtools' two BGZF threads mostly wait on it.
    bowtie2 -x "${index_name}" \\
        --threads ${task.cpus} \\
        ${aligner_args} \\
        -1 "${r1}" \\
        -2 "${r2}" \\
        2> "${prefix}.bowtie2.log" \\
        | samtools view -@ 2 -b -o "${prefix}.bam" -

    cat "${prefix}.bowtie2.log" >&2
    samtools quickcheck -v "${prefix}.bam"

    got=\$(awk '/ reads; of these:\$/ { print \$1; exit }' "${prefix}.bowtie2.log")
    if [ "\${got:-}" != "${n_pairs}" ]; then
        echo "error: bowtie2 read \${got:-no} pairs from batch ${batch} (L=${L}), but the" >&2
        echo "       batch plan says its FASTQ holds ${n_pairs}. The shard is incomplete." >&2
        exit 1
    fi
    echo "aligned ${n_pairs} pairs, matching the batch plan" >&2
    """

    stub:
    """
    touch "sim_L${L}_${batch}.bam" "sim_L${L}_${batch}.bowtie2.log"
    echo "bowtie2 stub" > bowtie2_version.txt
    """
}
