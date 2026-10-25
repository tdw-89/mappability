/*
 * Target loci: gene span + strand-aware flanks.
 *
 * Without --genes, bin/select_genes.sh picks every gene of --gene_biotype from
 * the GFF3. Biotypes are spelled differently by each annotation source
 * (gene_biotype= NCBI, biotype= Ensembl, gene_type= GENCODE), so it detects the
 * attribute, and when genes carry none it takes protein-coding genes to be the
 * ones with a CDS. gene_selection.tsv records which rule was applied. Genes on
 * sequences the FASTA lacks (alt haplotypes, patches) are left out and listed in
 * excluded_genes.tsv; if that leaves none, the run stops.
 *
 * IDs are taken from the GFF3 `ID=` attribute, not `Name=`: gene symbols repeat
 * (23 protein-coding symbols occur twice in zebrafish GRCz12tu), and
 * `unimap loci` skips an ambiguous alias rather than guessing, so selecting by
 * Name would silently drop those genes. With --genes the list is used as given.
 *
 * Every requested ID must come back as a locus. `unimap loci` warns about IDs
 * it cannot resolve and carries on; here a short loci.tsv stops the run, since
 * a gene that quietly drops out is a gene with no mappability track.
 *
 * Read names are `region|chrom|start|F|orient`, so a sequence name containing
 * `|` would only fail at scoring, after all the alignment; it is refused here.
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
    path "loci.tsv",           emit: loci
    path "gene_ids.txt",       emit: ids
    path "excluded_genes.tsv", emit: excluded
    path "gene_selection.tsv", emit: selection

    script:
    def project = params.unimap_project ? "export UNIMAP_PROJECT='${file(params.unimap_project)}'" : ''
    def select = genes.name == 'NO_GENES'
        ? "select_genes.sh \"${gff}\" \"${fai}\" '${params.gene_biotype}' '${params.biotype_attribute}' gene_ids.txt excluded_genes.tsv gene_selection.tsv"
        : "grep -v '^[[:space:]]*\$' \"${genes}\" > gene_ids.txt || true; " +
          "printf 'id\\tseqid\\treason\\n' > excluded_genes.tsv; " +
          "printf 'key\\tvalue\\nrule\\tlist %s\\n' '${genes.name}' > gene_selection.tsv"
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

    bad_seqs=\$(awk -F'\\t' 'NR > 1 && \$2 ~ /[|]/ && !seen[\$2]++ && n++ < 5 { print \$2 }' loci.tsv)
    if [ -n "\${bad_seqs}" ]; then
        echo "error: sequence names containing '|' cannot be used (UniMap read names are" >&2
        echo "       '|'-separated); rename them in both the FASTA and GFF:" >&2
        echo "\${bad_seqs}" >&2
        exit 1
    fi
    """

    stub:
    """
    printf 'locus_id\\n' > loci.tsv
    touch gene_ids.txt excluded_genes.tsv gene_selection.tsv
    """
}
