#!/usr/bin/env bash
#
# plan_batches.sh <regions.tsv> <tiling.tsv> <n_batches> <pairs_per_batch> <plan.tsv> <batches.tsv>
#
# Pack simulation regions into batches of near-equal read-pair count (greedy
# longest-first: biggest region to the lightest batch). Each batch becomes one
# ALIGN_BATCH task per read length, so the balance is what sets the wall time of
# the alignment step.
#
# How many batches: n_batches if it is > 0; otherwise just enough that the
# read length with the most pairs averages at most pairs_per_batch per batch,
# so task length stays the same whatever the genome and grid. Either way capped
# at 999 (batch IDs are three digits) and at the number of regions -- a region
# is never split, so one larger than pairs_per_batch makes a larger batch.
#
# Round-robin by region would balance the number of regions, not the work:
# regions run from a few kb to over a Mb, so a round-robin batch can carry
# several times the reads of another.
#
# One assignment serves every read length. A region's pair count differs
# between read lengths only through fragments skipped for Ns, so a balanced
# plan for the total is balanced for each L.
#
# Writes
#   plan.tsv     region_id, batch -- in regions.tsv order within each batch
#   batches.tsv  L, batch, n_regions, n_pairs -- the pairs each batch's FASTQ
#                must hold, which ALIGN_BATCH checks against bowtie2's count

set -euo pipefail
export LC_ALL=C

regions=$1 tiling=$2 n_batches=$3 pairs_per_batch=$4 plan=$5 manifest=$6

# region_id <tab> total pairs (all L) <tab> position in regions.tsv
awk -F'\t' '
    FNR == 1 { for (i = 1; i <= NF; i++) col[FILENAME, $i] = i; next }
    FILENAME == ARGV[1] { order[$col[FILENAME, "locus_id"]] = FNR - 1; next }
    { pairs[$col[FILENAME, "region_id"]] += $col[FILENAME, "n_pairs"] }
    END {
        for (r in order) {
            if (!(r in pairs)) { print "region " r " has no rows in tiling.tsv" > "/dev/stderr"; exit 1 }
            printf "%s\t%d\t%d\n", r, pairs[r], order[r]
        }
    }' "${regions}" "${tiling}" > region_sizes.tsv

n_regions=$(wc -l < region_sizes.tsv | tr -d ' ')
if [[ "${n_regions}" -eq 0 ]]; then
    echo "error: no regions in ${regions}" >&2
    exit 1
fi
if [[ "${n_batches}" -gt 0 ]]; then
    B=${n_batches}
else
    # The busiest read length's pairs, over the budget, rounded up.
    B=$(awk -F'\t' -v ppb="${pairs_per_batch}" '
        FNR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        { per_L[$col["L"]] += $col["n_pairs"] }
        END {
            m = 0; for (L in per_L) if (per_L[L] > m) m = per_L[L]
            b = int(m / ppb); if (b * ppb < m) b++
            print (b < 1 ? 1 : b)
        }' "${tiling}")
fi
B=$(( B < 999 ? B : 999 ))
B=$(( B < n_regions ? B : n_regions ))

# Largest first, ties by regions.tsv order, so the plan is deterministic.
sort -t $'\t' -k2,2nr -k3,3n region_sizes.tsv \
    | awk -F'\t' -v B="${B}" '
        BEGIN { for (b = 1; b <= B; b++) load[b] = 0 }
        {
            best = 1
            for (b = 2; b <= B; b++) if (load[b] < load[best]) best = b
            load[best] += $2
            printf "%s\t%03d\t%d\n", $1, best, $3
        }' \
    | sort -t $'\t' -k2,2 -k3,3n \
    | awk -F'\t' 'BEGIN { print "region_id\tbatch" } { print $1 "\t" $2 }' > "${plan}"

# Per (L, batch) pair counts, straight from tiling.tsv.
awk -F'\t' '
    FILENAME == ARGV[1] { if (FNR > 1) batch[$1] = $2; next }
    FNR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
    {
        r = $col["region_id"]; L = $col["L"]
        if (!(r in batch)) { print "region " r " in tiling.tsv but not in the plan" > "/dev/stderr"; exit 1 }
        key = L "\t" batch[r]
        pairs[key] += $col["n_pairs"]
        if (!((key, r) in seen)) { seen[key, r] = 1; nreg[key]++ }
    }
    END {
        print "L\tbatch\tn_regions\tn_pairs"
        for (k in pairs) print k "\t" nreg[k] "\t" pairs[k]
    }' "${plan}" "${tiling}" \
    | { IFS= read -r header; echo "${header}"; sort -t $'\t' -k1,1n -k2,2; } > "${manifest}"

rm -f region_sizes.tsv

echo "packed ${n_regions} region(s) into ${B} batch(es):" >&2
awk -F'\t' 'NR > 1 { printf "  L=%s batch %s: %d region(s), %d pair(s)\n", $1, $2, $3, $4 }' "${manifest}" >&2
