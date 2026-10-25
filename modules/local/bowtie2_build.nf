/*
 * Build the bowtie2 index when --bowtie2_index is not given.
 *
 * Meant for the test fixture and small genomes. A vertebrate genome takes hours
 * and tens of GB here; build that once and pass --bowtie2_index instead.
 */

process BOWTIE2_BUILD {
    tag   { fasta.name }

    container params.container_bowtie2

    input:
    path fasta
    val  index_name

    output:
    path "${index_name}.*.bt2*", emit: index

    script:
    """
    bowtie2-build --threads ${task.cpus} --seed 1 "${fasta}" "${index_name}"
    """

    stub:
    """
    for s in 1 2 3 4 rev.1 rev.2; do touch "${index_name}.\${s}.bt2"; done
    """
}
