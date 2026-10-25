#!/usr/bin/env bash
#
# select_genes.sh <annotation.gff3[.gz]> <ref.fa.fai> <biotypes> <attribute> \
#                 <gene_ids.txt> <excluded_genes.tsv> <gene_selection.tsv>
#
# Pick the genes of the given biotype(s) from a GFF3, by ID=, for `unimap loci`.
#
# GFF3 has no standard biotype attribute, so the rule depends on the source:
#   NCBI RefSeq    gene_biotype=protein_coding
#   Ensembl        biotype=protein_coding         (ID=gene:ENSG...)
#   GENCODE        gene_type=protein_coding
#   others         often none (MAKER, BRAKER, FlyBase-style): a gene is
#                  protein-coding if a CDS descends from it
# <attribute> "auto" takes whichever of gene_biotype, biotype, gene_type the
# gene features carry (the most common, if several). If no gene carries any of
# them, genes are selected by CDS descent instead -- which can only answer
# "protein_coding". Any other attribute name is used as given and must exist.
#
# A gene is a root feature (an ID and no Parent) whose type contains "gene",
# the same features `unimap loci` treats as genes. CDS descent follows Parent=
# up to the root; CDS marked pseudo=true and roots whose type contains "pseudo"
# do not count.
#
# Genes on sequences absent from the FASTA (alt haplotypes, patches, unplaced
# scaffolds a primary-assembly FASTA leaves out) cannot be simulated. They are
# written to <excluded_genes.tsv> and warned about rather than failing the run,
# unless none remain -- which means the GFF and FASTA name sequences
# differently (chr1 vs 1 vs NC_...).
#
# Writes
#   gene_ids.txt        selected IDs, in GFF order
#   excluded_genes.tsv  id, seqid, reason
#   gene_selection.tsv  key/value summary of the rule that was applied

set -euo pipefail
export LC_ALL=C

gff=$1 fai=$2 biotypes=$3 attribute=$4 out_ids=$5 out_excluded=$6 out_summary=$7

gzip -cdf "${gff}" | awk -F'\t' \
    -v biotypes="${biotypes}" -v attribute="${attribute}" \
    -v out_ids="${out_ids}" -v out_excluded="${out_excluded}" -v out_summary="${out_summary}" '
    # Fill A[key] = value from a GFF3 attribute column.
    function parse(attrs,   n, i, f, eq) {
        delete A
        n = split(attrs, f, ";")
        for (i = 1; i <= n; i++) {
            sub(/^[ \t]+/, "", f[i])
            eq = index(f[i], "=")
            if (eq > 1) A[substr(f[i], 1, eq - 1)] = substr(f[i], eq + 1)
        }
    }
    function fail(msg) { print "error: " msg > "/dev/stderr"; failed = 1; exit 1 }

    BEGIN {
        nb = split(biotypes, bt, ",")
        for (i = 1; i <= nb; i++) { gsub(/^[ \t]+|[ \t]+$/, "", bt[i]); want[bt[i]] = 1 }
        ncand = split("gene_biotype,biotype,gene_type", cand, ",")
        if (attribute != "auto") { ncand = 1; cand[1] = attribute }
    }

    # The .fai: sequences the FASTA has.
    FNR == NR { inseq[$1] = 1; next }

    /^##FASTA/ { exit }          # an embedded FASTA section ends the features
    /^#/       { next }
    { sub(/\r$/, "") }
    NF < 9     { next }
    {
        parse($9)
        type = $3
        id = ("ID" in A) ? A["ID"] : ""
        par = ("Parent" in A) ? A["Parent"] : ""
        if (type == "CDS") {
            if (par != "" && A["pseudo"] != "true") {
                n = split(par, ps, ",")
                for (i = 1; i <= n; i++) cdspar[ps[i]] = 1
            }
            next
        }
        if (type == "exon" || id == "") next
        if (par != "") { split(par, ps, ","); parent[id] = ps[1]; next }
        if (id in rtype) next    # a feature split over several lines
        rtype[id] = type; rseq[id] = $1; order[++nroot] = id
        if (tolower(type) !~ /gene/) next
        ngene++
        for (k = 1; k <= ncand; k++) {
            if (cand[k] in A) { bio[id, cand[k]] = A[cand[k]]; keycount[cand[k]]++ }
        }
    }

    END {
        if (failed) exit 1
        if (ngene == 0) fail("no gene features (a root feature whose type contains \"gene\") in the annotation")

        key = ""
        for (k = 1; k <= ncand; k++)
            if (keycount[cand[k]] > 0 && (key == "" || keycount[cand[k]] > keycount[key])) key = cand[k]
        if (attribute != "auto" && key == "")
            fail("no gene feature has a " attribute "= attribute")
        if (key != "") {
            rule = "biotype"
        } else {
            rule = "cds"
            if (!(nb == 1 && bt[1] == "protein_coding"))
                fail("genes carry no gene_biotype, biotype or gene_type attribute, so only " \
                     "protein_coding can be selected (by CDS descent); got --gene_biotype " biotypes \
                     ". Name the attribute with --biotype_attribute if it is called something else.")
            for (p in cdspar) {
                r = p
                for (i = 0; i < 100 && (r in parent); i++) r = parent[r]
                if (r in rtype) cdsroot[r] = 1
            }
        }

        printf "" > out_ids
        printf "id\tseqid\treason\n" > out_excluded
        for (j = 1; j <= nroot; j++) {
            id = order[j]; t = tolower(rtype[id])
            if (t !~ /gene/) continue
            if (rule == "biotype") {
                if (!((id, key) in bio)) continue
                seen[bio[id, key]]++
                if (!(bio[id, key] in want)) continue
            } else if (!(id in cdsroot) || t ~ /pseudo/) {
                continue
            }
            if (!(rseq[id] in inseq)) {
                printf "%s\t%s\tsequence not in FASTA\n", id, rseq[id] > out_excluded
                nexcl++; exseq[rseq[id]] = 1
                continue
            }
            print id > out_ids
            nsel++
        }

        printf "key\tvalue\n" > out_summary
        printf "rule\t%s\n", (rule == "biotype" ? "attribute " key : "CDS descent") > out_summary
        printf "biotypes\t%s\n", biotypes > out_summary
        printf "gene_features\t%d\n", ngene > out_summary
        printf "selected\t%d\n", nsel > out_summary
        printf "excluded_no_sequence\t%d\n", nexcl > out_summary

        if (rule == "biotype")
            printf "selecting genes with %s in {%s} (%d gene features)\n", key, biotypes, ngene > "/dev/stderr"
        else
            printf "no biotype attribute on genes: selecting genes with a CDS (%d gene features)\n", ngene > "/dev/stderr"
        if (nexcl > 0) {
            ns = 0; for (s in exseq) ns++
            printf "warning: %d selected gene(s) on %d sequence(s) not in the FASTA were left out; see %s\n", \
                   nexcl, ns, out_excluded > "/dev/stderr"
        }
        if (nsel == 0) {
            if (nexcl > 0)
                fail("every selected gene is on a sequence missing from the FASTA: the GFF and FASTA " \
                     "name sequences differently (e.g. chr1 vs 1)")
            msg = "no gene matched --gene_biotype " biotypes
            if (rule == "biotype") {
                msg = msg "; values of " key " present:"
                for (v in seen) msg = msg "\n  " v " (" seen[v] ")"
            }
            fail(msg)
        }
        printf "selected %d gene(s)\n", nsel > "/dev/stderr"
    }' "${fai}" -
