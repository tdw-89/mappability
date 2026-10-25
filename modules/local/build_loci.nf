/*
 * Target loci: gene span + strand-aware flanks.
 *
 * IDs are taken from the GFF3 `ID=` attribute, not `Name=`: in GRCz12tu 23
 * protein-coding gene symbols occur more than once, and `unimap loci` skips an
 * ambiguous alias rather than guessing, so selecting by Name would silently
 * drop those genes. With --genes the list is used as given.
 *
 * Every requested ID must come back as a locus. `unimap loci` warns about IDs
 * it cannot resolve and carries on; here a short loci.tsv stops the run, since
 * a gene that quietly drops out is a gene with no mappability track.
 */

process BUILD_LOCI {
    tag   { params.genes ? file(params.genes).name : params.gene_biotype }

    container params.container_unimap

    publishDir "${params.outdir}/loci", mode: params.publish_mode

    input:
    path fasta
    path fai
    path gff
    path genes      // assets/NO_GENES unless --genes was given

    output:
    path "loci.tsv",   emit: loci
    path "gene_ids.txt", emit: ids

    script:
    def project = params.unimap_project ? "export UNIMAP_PROJECT='${file(params.unimap_project)}'" : ''
    def select = genes.name == 'NO_GENES'
        ? """\
        gzip -cdf "${gff}" | awk -F'\\t' -v biotype='${params.gene_biotype}' '
            \$1 !~ /^#/ && \$3 == "gene" && \$9 ~ ("(^|;)gene_biotype=" biotype "(;|\$)") {
                n = split(\$9, attrs, ";")
                for (i = 1; i <= n; i++) {
                    if (attrs[i] ~ /^ID=/) {
                        print substr(attrs[i], 4)
                        break
                    }
                }
            }' > gene_ids.txt
        """.stripIndent()
        : """\
        grep -v '^[[:space:]]*\$' "${genes}" > gene_ids.txt || true
        """.stripIndent()
    """
    ${project}
    export UNIMAP_THREADS=${params.julia_threads ?: task.cpus}

    ${select}
    n_ids=\$(sort -u gene_ids.txt | wc -l | tr -d ' ')
    if [ "\${n_ids}" -eq 0 ]; then
        echo "error: no gene IDs selected" >&2
        exit 1
    fi
    echo "selected \${n_ids} gene ID(s)" >&2

    run_unimap loci \\
        -f "${fasta}" \\
        -a "${gff}" \\
        -g gene_ids.txt \\
        --upstream ${params.upstream} \\
        --downstream ${params.downstream} \\
        -o loci.tsv

    n_loci=\$(( \$(wc -l < loci.tsv) - 1 ))
    if [ "\${n_loci}" -ne "\${n_ids}" ]; then
        echo "error: \${n_ids} gene ID(s) requested but only \${n_loci} locus/loci built;" >&2
        echo "       the IDs unimap could not resolve are in its warnings above." >&2
        exit 1
    fi
    echo "built \${n_loci} loci" >&2
    """

    stub:
    """
    printf 'locus_id\\n' > loci.tsv
    touch gene_ids.txt
    """
}
