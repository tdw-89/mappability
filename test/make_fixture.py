#!/usr/bin/env python3
"""
Build the test fixture: a ~120 kb three-chromosome genome and a GFF3 laid out
so each feature of the pipeline has something to bite on.

    test/make_fixture.py test/data        # then: bowtie2-build (see run_tests.sh)

Deterministic (fixed seed), standard library only. The layout:

  gene-dupA     chr1, 3 kb of which is copied verbatim onto chr2. Fragments
                from the copy cannot be placed uniquely: low MAPQ, ex-mapping.
  gene-uniqB    chr1, unique sequence: map-back ~1 at every cutoff.
  gene-nearC    chr1, 1 kb downstream of uniqB, so the two loci merge into one
                simulation region (gap <= max F - 2).
  gene-edgeD    chr1, starts 200 bp from the chromosome start: flank clipping.
  gene-revE     chr2, minus strand.
  gene-nblkF    chr3, with a 500 bp run of N inside it: N-skipped fragments.
  gene-sameG    chr3, shares Name= with uniqB, so selecting by Name would be
                ambiguous; selecting by ID (what BUILD_LOCI does) is not.
  rna-lnc1      chr2, lncRNA: must NOT be selected.
  gene-pseudo1  chr1, pseudogene: must NOT be selected.
"""

import os
import random
import sys

SEED = 20261024
LINE = 80


def rand_seq(rng, n):
    return [rng.choice("ACGT") for _ in range(n)]


def main(outdir):
    rng = random.Random(SEED)
    os.makedirs(outdir, exist_ok=True)

    chr1 = rand_seq(rng, 60_000)
    chr2 = rand_seq(rng, 40_000)
    chr3 = rand_seq(rng, 20_000)

    # Paralog: chr1[10_000:13_000) copied to chr2[20_000:23_000).
    chr2[20_000:23_000] = chr1[10_000:13_000]
    # Assembly gap inside gene-nblkF.
    chr3[5_000:5_500] = ["N"] * 500

    with open(os.path.join(outdir, "genome.fa"), "w") as fa:
        for name, seq in (("chr1", chr1), ("chr2", chr2), ("chr3", chr3)):
            fa.write(f">{name} test fixture\n")
            s = "".join(seq)
            for i in range(0, len(s), LINE):
                fa.write(s[i:i + LINE] + "\n")

    # (id, chrom, start1, end1, strand, biotype, name, exons as (start1, end1))
    genes = [
        ("dupA",    "chr1", 10_500, 12_500, "+", "protein_coding", "dupa",    [(10_500, 11_000), (12_000, 12_500)]),
        ("uniqB",   "chr1", 30_001, 33_000, "-", "protein_coding", "shared",  [(30_001, 31_000), (32_000, 33_000)]),
        ("nearC",   "chr1", 34_001, 35_000, "+", "protein_coding", "nearc",   [(34_001, 35_000)]),
        ("edgeD",   "chr1",    201,  1_200, "+", "protein_coding", "edged",   [(201, 1_200)]),
        ("pseudo1", "chr1", 45_001, 46_000, "+", "pseudogene",     "pseudo1", [(45_001, 46_000)]),
        ("revE",    "chr2",  5_001,  7_000, "-", "protein_coding", "reve",    [(5_001, 5_600), (6_400, 7_000)]),
        ("lnc1",    "chr2", 30_001, 31_000, "+", "lncRNA",         "lnc1",    [(30_001, 31_000)]),
        ("nblkF",   "chr3",  4_501,  6_500, "+", "protein_coding", "nblkf",   [(4_501, 6_500)]),
        ("sameG",   "chr3", 15_001, 17_000, "-", "protein_coding", "shared",  [(15_001, 17_000)]),
    ]

    with open(os.path.join(outdir, "annotation.gff3"), "w") as gff:
        gff.write("##gff-version 3\n")
        for name, seq in (("chr1", chr1), ("chr2", chr2), ("chr3", chr3)):
            gff.write(f"##sequence-region {name} 1 {len(seq)}\n")
        for gid, chrom, s, e, strand, biotype, name, exons in genes:
            prefix = "rna" if biotype == "lncRNA" else "gene"
            tx_type = "lnc_RNA" if biotype == "lncRNA" else "mRNA"
            gene_id = f"{prefix}-{gid}"
            gff.write(f"{chrom}\tTest\tgene\t{s}\t{e}\t.\t{strand}\t.\t"
                      f"ID={gene_id};Name={name};gene_biotype={biotype}\n")
            if biotype == "pseudogene":
                continue
            tx = f"rna-{gid}.1"
            gff.write(f"{chrom}\tTest\t{tx_type}\t{s}\t{e}\t.\t{strand}\t.\t"
                      f"ID={tx};Parent={gene_id};Name={gid}.1;gene={name}\n")
            for k, (es, ee) in enumerate(exons, 1):
                gff.write(f"{chrom}\tTest\texon\t{es}\t{ee}\t.\t{strand}\t.\t"
                          f"ID=exon-{gid}.1-{k};Parent={tx}\n")

    with open(os.path.join(outdir, "genes_subset.txt"), "w") as f:
        f.write("gene-uniqB\ngene-dupA\n")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "data"))
