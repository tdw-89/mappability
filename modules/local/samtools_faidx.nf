/*
 * Index the reference when no .fai sits beside it. UniMap reads sequence
 * through the .fai and refuses a FASTA without one.
 */

process SAMTOOLS_FAIDX {
    tag   { fasta.name }

    container params.container_unimap

    input:
    path fasta

    output:
    path "${fasta}.fai", emit: fai

    script:
    """
    samtools faidx "${fasta}"
    """

    stub:
    """
    touch "${fasta}.fai"
    """
}
