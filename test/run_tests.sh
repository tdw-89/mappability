#!/bin/bash
#
# Tests for the mappability pipeline.
#
#   test/run_tests.sh              unit tests + launch checks (seconds; needs
#                                  nextflow, julia)
#   test/run_tests.sh --e2e        also the full pipeline on the fixture, compared
#                                  exactly against UniMap's own monolithic path
#                                  (needs bowtie2, samtools; ~3 min)
#   test/run_tests.sh --docker     also -profile test,docker (needs the
#                                  container images; slow under emulation)
#
# UniMap comes from $UNIMAP_PROJECT, default ../UniMap.jl beside this repo.
#
# The e2e comparison is the test that matters most: it runs the same loci and
# grid through `unimap simulate -> align-script -> align.sh -> score` (one
# bowtie2 run per read length, as UniMap intends) and requires every Zarr array
# and Parquet table to be identical to the pipeline's sharded, concatenated
# result. That is the claim the whole design rests on.

set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${HERE}/.." && pwd)"
BIN="${ROOT}/bin"
UNIMAP="${UNIMAP_PROJECT:-$(cd "${ROOT}/../UniMap.jl" 2>/dev/null && pwd)}"
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT

E2E=0
DOCKER=0
for arg in "$@"; do
    case "${arg}" in
        --e2e)    E2E=1 ;;
        --docker) E2E=1; DOCKER=1 ;;
        *) echo "unknown argument: ${arg}" >&2; exit 2 ;;
    esac
done

PASS=0
FAIL=0
FAILURES=()

ok()   { printf '  \033[32mok\033[0m   %s\n' "$1"; PASS=$((PASS + 1)); }
bad()  { printf '  \033[31mFAIL\033[0m %s\n' "$1"; FAIL=$((FAIL + 1)); FAILURES+=("$1"); }
assert_eq() { if [[ "$2" == "$3" ]]; then ok "$1"; else bad "$1 -- expected '$2', got '$3'"; fi; }
assert_contains() { if [[ "$3" == *"$2"* ]]; then ok "$1"; else bad "$1 -- '$2' not in output"; fi; }

# Output of test/check_scores.jl ("  ok   ..." / "  FAIL ..."), counted as tests.
tally() {
    local out="$1" n_ok n_fail
    printf '%s\n' "${out}"
    n_ok=$(printf '%s\n' "${out}" | grep -c '^  ok ')
    n_fail=$(printf '%s\n' "${out}" | grep -c '^  FAIL ')
    PASS=$((PASS + n_ok)); FAIL=$((FAIL + n_fail))
    if [[ ${n_fail} -gt 0 ]]; then
        while IFS= read -r line; do FAILURES+=("${line#  FAIL }"); done < <(printf '%s\n' "${out}" | grep '^  FAIL ')
    fi
    if [[ $((n_ok + n_fail)) -eq 0 ]]; then
        bad "check_scores.jl produced no results"
    fi
}

if [[ -z "${UNIMAP}" || ! -f "${UNIMAP}/Project.toml" ]]; then
    echo "error: set UNIMAP_PROJECT to a UniMap.jl checkout" >&2
    exit 2
fi
# Nextflow 25.04 runs on Java 17-24; a Mac may default to something newer.
if [[ -z "${JAVA_HOME:-}" && -x /usr/libexec/java_home ]]; then
    for v in 24 21 17; do
        if j="$(/usr/libexec/java_home -v "${v}" 2>/dev/null)"; then export JAVA_HOME="${j}"; break; fi
    done
fi

echo "mappability tests"
echo

# ---------------------------------------------------------------------------
echo "plan_batches.sh -- packing regions into balanced batches"

# Five regions, two read lengths. Pairs per region per L: r1 100, r2 60, r3 50,
# r4 30, r5 10 (r5 loses 2 pairs at L=50 to Ns). Totals across L: 200, 120,
# 100, 60, 18. Longest-first into 2 batches:
#   r1 -> b1 (200)   r2 -> b2 (120)   r3 -> b2 (220)   r4 -> b1 (260)   r5 -> b2 (238)
mkdir -p "${WORK}/plan"
cd "${WORK}/plan"
printf 'locus_id\tchrom\tstart0\n' > regions.tsv
for r in r1 r2 r3 r4 r5; do printf '%s\tchr1\t0\n' "${r}" >> regions.tsv; done
{
    printf 'region_id\tL\tF\tstep\tfirst_start0\tlast_start0\tn_starts\tn_skipped_N\tn_pairs\n'
    for spec in r1:100 r2:60 r3:50 r4:30 r5:10; do
        r=${spec%%:*}; n=${spec##*:}
        # two F rows per (region, L), splitting the pairs
        for L in 30 50; do
            m=${n}; [[ "${r}" == r5 && ${L} == 50 ]] && m=$((n - 2))
            printf '%s\t%s\t80\t5\t0\t0\t0\t0\t%s\n'  "${r}" "${L}" $((m / 2))
            printf '%s\t%s\t200\t5\t0\t0\t0\t0\t%s\n' "${r}" "${L}" $((m - m / 2))
        done
    done
} > tiling.tsv

"${BIN}/plan_batches.sh" regions.tsv tiling.tsv 2 0 plan.tsv batches.tsv 2> /dev/null
assert_eq "planner exits cleanly"              "0" "$?"
assert_eq "every region planned exactly once"  "r1 r2 r3 r4 r5" "$(tail -n +2 plan.tsv | cut -f1 | sort | tr '\n' ' ' | sed 's/ $//')"
assert_eq "longest-first assignment"           "r1:001 r2:002 r3:002 r4:001 r5:002" \
    "$(tail -n +2 plan.tsv | sort | awk -F'\t' '{printf "%s%s:%s", (NR>1?" ":""), $1, $2}')"
assert_eq "regions keep regions.tsv order within a batch" "r2 r3 r5" \
    "$(awk -F'\t' '$2=="002"{print $1}' plan.tsv | tr '\n' ' ' | sed 's/ $//')"
assert_eq "manifest rows: 2 L x 2 batches"     "4" "$(tail -n +2 batches.tsv | wc -l | tr -d ' ')"
assert_eq "L=30 batch 001 = r1 + r4 = 130 pairs" "130" "$(awk -F'\t' '$1==30 && $2=="001"{print $4}' batches.tsv)"
assert_eq "L=50 batch 002 = r2 + r3 + r5 = 118"  "118" "$(awk -F'\t' '$1==50 && $2=="002"{print $4}' batches.tsv)"
for L in 30 50; do
    assert_eq "L=${L}: manifest pairs sum to tiling.tsv" \
        "$(awk -F'\t' -v L=${L} 'NR>1 && $2==L {s+=$9} END{print s}' tiling.tsv)" \
        "$(awk -F'\t' -v L=${L} 'NR>1 && $1==L {s+=$4} END{print s}' batches.tsv)"
done

"${BIN}/plan_batches.sh" regions.tsv tiling.tsv 40 0 plan40.tsv batches40.tsv 2> /dev/null
assert_eq "n_batches capped at the region count" "5" "$(tail -n +2 plan40.tsv | cut -f2 | sort -u | wc -l | tr -d ' ')"

# By pair budget: the busiest L has 250 pairs.
n_planned() { tail -n +2 "$1" | cut -f2 | sort -u | wc -l | tr -d ' '; }
"${BIN}/plan_batches.sh" regions.tsv tiling.tsv 0 100 pa.tsv ba.tsv 2> /dev/null
assert_eq "pairs_per_batch 100: ceil(250 / 100) = 3 batches" "3" "$(n_planned pa.tsv)"
"${BIN}/plan_batches.sh" regions.tsv tiling.tsv 0 250 pa.tsv ba.tsv 2> /dev/null
assert_eq "pairs_per_batch 250: exactly one batch's worth = 1" "1" "$(n_planned pa.tsv)"
"${BIN}/plan_batches.sh" regions.tsv tiling.tsv 0 1 pa.tsv ba.tsv 2> /dev/null
assert_eq "pairs_per_batch 1: capped at the region count" "5" "$(n_planned pa.tsv)"

grep -v '^r3' tiling.tsv > tiling_missing.tsv
"${BIN}/plan_batches.sh" regions.tsv tiling_missing.tsv 2 0 p.tsv b.tsv 2> err.txt
rc=$?
if [[ ${rc} -ne 0 ]]; then ok "region absent from tiling.tsv is an error"; else bad "region absent from tiling.tsv is an error"; fi
cd "${ROOT}"

# ---------------------------------------------------------------------------
echo
echo "select_genes.sh -- protein-coding genes from any GFF3"

# The fixture's GFF is NCBI-style (gene_biotype=, ID=gene-...). Rewrite it in the
# other conventions; each must give the same 7 protein-coding genes and leave
# out the lncRNA and the pseudogene.
G="${WORK}/genes"
mkdir -p "${G}"
cd "${G}"
fixture_gff="${ROOT}/test/data/annotation.gff3"
# A .fai (names and lengths are all select_genes.sh and unimap loci need of it
# here, but make it a real one).
cp "${ROOT}/test/data/genome.fa" genome.fa
awk 'function out() { if (name != "") printf "%s\t%d\t%d\t%d\t%d\n", name, len, off, lb, lb + 1 }
     /^>/ { out(); name = substr($1, 2); len = 0; lb = 0; pos += length($0) + 1; off = pos; next }
     { if (lb == 0) lb = length($0); len += length($0); pos += length($0) + 1 }
     END { out() }' genome.fa > genome.fa.fai

cp "${fixture_gff}" ncbi.gff3
# An embedded FASTA section ends the features: a gene after it is not one.
printf '##FASTA\nchr1\tX\tgene\t1\t10\t.\t+\t.\tID=gene-bogus;gene_biotype=protein_coding\n' >> ncbi.gff3
sed -e 's/ID=gene-/ID=gene:/; s/Parent=gene-/Parent=gene:/; s/gene_biotype=/biotype=/' "${fixture_gff}" > ensembl.gff3
sed -e 's/gene_biotype=/gene_type=/' "${fixture_gff}" | gzip -c > gencode.gff3.gz
# No biotype attribute at all: protein-coding genes are the ones with a CDS.
# The pseudogene gets a CDS under a "pseudogene"-typed gene, the lncRNA a
# pseudo=true CDS; neither may count. CRLF line endings, as some tools write.
awk -F'\t' -v OFS='\t' '
    $3 == "gene" && $9 ~ /gene_biotype=protein_coding/ { pc[substr($9, 4, index($9, ";") - 4)] = 1 }
    $3 == "mRNA" { split($9, a, ";"); id = substr(a[1], 4); par = substr(a[2], 8); if (par in pc) cds[id] = 1 }
    { line = $0; sub(/;gene_biotype=[^;]*/, "", line) }
    $9 ~ /^ID=gene-pseudo1/ { sub(/\tgene\t/, "\tpseudogene\t", line) }
    { print line }
    $3 == "exon" { par = $9; sub(/.*Parent=/, "", par); sub(/;.*/, "", par)
                   if (par in cds) { $3 = "CDS"; $9 = "Parent=" par; print } }
    $9 ~ /^ID=gene-pseudo1/ { $3 = "CDS"; $9 = "ID=cds-p1;Parent=gene-pseudo1"; print }
    $9 ~ /^ID=rna-lnc1\.1/  { $3 = "CDS"; $9 = "ID=cds-l1;Parent=rna-lnc1.1;pseudo=true"; print }
' "${fixture_gff}" | awk '{ printf "%s\r\n", $0 }' > nobiotype.gff3

expected_ids="dupA edgeD nblkF nearC revE sameG uniqB"
ids_of() { sed -e 's/^gene[-:]//' "$1" | sort | tr '\n' ' ' | sed 's/ $//'; }
for style in ncbi:ncbi.gff3:gene_biotype ensembl:ensembl.gff3:biotype gencode:gencode.gff3.gz:gene_type \
             nobiotype:nobiotype.gff3:CDS; do
    IFS=: read -r name gff rule <<< "${style}"
    "${BIN}/select_genes.sh" "${gff}" genome.fa.fai protein_coding auto \
        "${name}_ids.txt" "${name}_excl.tsv" "${name}_sel.tsv" 2> "${name}.err"
    assert_eq "${name}: exits cleanly" "0" "$?"
    assert_eq "${name}: the 7 protein-coding genes, no lncRNA/pseudogene" "${expected_ids}" "$(ids_of "${name}_ids.txt")"
    assert_contains "${name}: rule recorded" "${rule}" "$(grep '^rule' "${name}_sel.tsv")"
done
assert_eq "ensembl: IDs kept as written (gene:...)" "gene:dupA" "$(head -n 1 ensembl_ids.txt)"

"${BIN}/select_genes.sh" ncbi.gff3 genome.fa.fai protein_coding,lncRNA auto l.txt e.tsv s.tsv 2> /dev/null
assert_eq "biotype list: protein_coding,lncRNA = 8" "8" "$(wc -l < l.txt | tr -d ' ')"

grep -v '^chr3' genome.fa.fai > no_chr3.fai
"${BIN}/select_genes.sh" ncbi.gff3 no_chr3.fai protein_coding auto m.txt m_excl.tsv m_sel.tsv 2> m.err
assert_eq "genes on sequences not in the FASTA are left out" "dupA edgeD nearC revE uniqB" "$(ids_of m.txt)"
assert_eq "and listed in excluded_genes.tsv" "gene-nblkF gene-sameG" "$(tail -n +2 m_excl.tsv | cut -f1 | sort | tr '\n' ' ' | sed 's/ $//')"
assert_contains "and warned about" "2 selected gene(s) on 1 sequence(s) not in the FASTA" "$(cat m.err)"

sed 's/^chr//' genome.fa.fai > renamed.fai
out=$("${BIN}/select_genes.sh" ncbi.gff3 renamed.fai protein_coding auto x.txt x.tsv xs.tsv 2>&1); rc=$?
assert_eq "GFF and FASTA naming sequences differently is an error" "1" "${rc}"
assert_contains "that says so" "name sequences differently" "${out}"
out=$("${BIN}/select_genes.sh" ncbi.gff3 genome.fa.fai protein-coding auto x.txt x.tsv xs.tsv 2>&1); rc=$?
assert_eq "a biotype no gene has is an error" "1" "${rc}"
assert_contains "that lists the values present" "lncRNA (1)" "${out}"
out=$("${BIN}/select_genes.sh" nobiotype.gff3 genome.fa.fai lncRNA auto x.txt x.tsv xs.tsv 2>&1); rc=$?
assert_eq "without a biotype attribute, only protein_coding can be selected" "1" "${rc}"
out=$("${BIN}/select_genes.sh" ncbi.gff3 genome.fa.fai protein_coding gene_kind x.txt x.tsv xs.tsv 2>&1); rc=$?
assert_eq "a named --biotype_attribute that no gene has is an error" "1" "${rc}"

# The selected IDs must resolve in `unimap loci` -- the run stops otherwise.
for name in ensembl nobiotype; do
    gff=$(ls "${name}".gff3*)
    julia --startup-file=no --project="${UNIMAP}" -m UniMap loci -f genome.fa -a "${gff}" \
        -g "${name}_ids.txt" --upstream 500 --downstream 500 -o "${name}_loci.tsv" > "${name}_loci.log" 2>&1
    assert_eq "${name}: unimap loci builds all 7 loci from those IDs" "7" "$(( $(wc -l < "${name}_loci.tsv" 2>/dev/null || echo 1) - 1 ))"
done
cd "${ROOT}"

# ---------------------------------------------------------------------------
echo
echo "check_score.jl -- every simulated pair scored once"

cat > "${WORK}/simulate.toml" <<'EOF'
[output.n_pairs_per_L]
30 = 1000
50 = 990
EOF
cat > "${WORK}/score_ok.toml" <<'EOF'
[bam_stats.30]
pairs = 1000
incomplete_pairs = 0
records = 2000
[bam_stats.50]
pairs = 990
incomplete_pairs = 0
records = 1980
EOF
sed 's/^pairs = 990/pairs = 900/' "${WORK}/score_ok.toml" > "${WORK}/score_short.toml"
awk '!done && /incomplete_pairs = 0/ { sub(/= 0/, "= 3"); done = 1 } 1' \
    "${WORK}/score_ok.toml" > "${WORK}/score_incomplete.toml"

"${BIN}/check_score.jl" "${WORK}/simulate.toml" "${WORK}/score_ok.toml" 2> /dev/null
assert_eq "complete BAMs pass" "0" "$?"
out=$("${BIN}/check_score.jl" "${WORK}/simulate.toml" "${WORK}/score_short.toml" 2>&1); rc=$?
assert_eq "a missing shard's pairs fail the run" "1" "${rc}"
assert_contains "and say how many were lost" "simulated 990 pairs, scored 900" "${out}"
"${BIN}/check_score.jl" "${WORK}/simulate.toml" "${WORK}/score_incomplete.toml" 2> /dev/null
assert_eq "pairs with a missing mate fail the run" "1" "$?"

# ---------------------------------------------------------------------------
echo
echo "launch-time validation"

nf_launch() {   # runs from WORK so no .nextflow/ or results land in the repo
    (cd "${WORK}" && nextflow run "${ROOT}" -profile test,local --unimap_project "${UNIMAP}" \
        -ansi-log false "$@" 2>&1)
}

out=$(nf_launch --bowtie2_args '--very-sensitive --no-unal' -stub-run)
assert_contains "--no-unal in --bowtie2_args is refused" "must not set --no-unal" "${out}"
out=$(nf_launch --bowtie2_args '-X 300' -stub-run)
assert_contains "-X in --bowtie2_args is refused (it is max F)" "must not set -X" "${out}"
# Without the test profile there is no grid: it must be given.
out=$(cd "${WORK}" && nextflow run "${ROOT}" -profile local --unimap_project "${UNIMAP}" -ansi-log false -stub-run \
          --fasta "${ROOT}/test/data/genome.fa" --gff "${ROOT}/test/data/annotation.gff3" 2>&1)
assert_contains "the grid is required" "missing required parameter(s): --read_lengths, --fragment_lengths" "${out}"
out=$(nf_launch --read_lengths -stub-run)
assert_contains "a bare --read_lengths is missing, not 'true'" "missing required parameter(s): --read_lengths" "${out}"
out=$(nf_launch --fragment_lengths '50,abc' -stub-run)
assert_contains "malformed --fragment_lengths is refused" "comma-separated list" "${out}"
out=$(nf_launch --bowtie2_index /nonexistent/idx -stub-run)
assert_contains "incomplete --bowtie2_index is refused" "no complete bowtie2 index" "${out}"
out=$(cd "${WORK}" && nextflow run "${ROOT}" -profile test -ansi-log false -stub-run 2>&1)
assert_contains "no container and no --unimap_project is refused" "--unimap_project" "${out}"
# A cloud batch executor runs every task in its container even with no engine
# enabled -- how Seqera's AWS Batch environments are configured.
printf 'process.executor = "awsbatch"\n' > "${WORK}/awsbatch.config"
out=$(cd "${WORK}" && nextflow run "${ROOT}" -profile test -c "${WORK}/awsbatch.config" -ansi-log false -preview 2>&1)
if [[ "${out}" != *"--unimap_project"* ]]; then ok "awsbatch with no container engine is accepted"; else bad "awsbatch with no container engine is accepted"; printf '%s\n' "${out}" | tail -5; fi
# Under Fusion the work dir is S3; the two tasks that write throwaway files
# must run in local scratch, without losing the profile's resource blocks.
out=$(nextflow config "${ROOT}" -profile seqera -flat 2>&1)
assert_contains "seqera: SIMULATE_AND_BATCH and SCORE use local scratch" \
    "process.'withName:SIMULATE_AND_BATCH|SCORE'.scratch = true" "${out}"
assert_contains "seqera: and keep their resources" \
    "process.'withName:SIMULATE_AND_BATCH'.cpus = 32" "${out}"
out=$(nf_launch -stub-run --outdir "${WORK}/stub")
assert_contains "stub run wires every process" "> SCORE" "${out}"
if [[ "${out}" == *"Done in"* ]]; then ok "stub run completes"; else bad "stub run completes"; printf '%s\n' "${out}" | tail -20; fi

# ---------------------------------------------------------------------------
if [[ ${E2E} -eq 1 ]]; then
    echo
    echo "end to end -- sharded pipeline vs monolithic UniMap (test,local)"

    for tool in bowtie2 bowtie2-build samtools; do
        command -v "${tool}" > /dev/null || { echo "error: ${tool} not on PATH" >&2; exit 2; }
    done

    # One index for both paths, so any difference is the sharding.
    mkdir -p "${WORK}/index"
    bowtie2-build -q --seed 1 "${ROOT}/test/data/genome.fa" "${WORK}/index/genome" > /dev/null

    # With the scratch settings the seqera profile uses, so the unstaging they
    # rely on (only declared outputs leave the task) is checked end to end.
    out=$(nf_launch --bowtie2_index "${WORK}/index/genome" --outdir "${WORK}/nf" -w "${WORK}/work" \
              -c "${ROOT}/conf/fusion_scratch.config")
    if [[ "${out}" == *"Done in"* ]]; then
        ok "pipeline completes"
    else
        bad "pipeline completes"; printf '%s\n' "${out}" | tail -30
    fi
    # Each of the two logs the dir it ran in; neither may be under the work dir.
    task_dirs=$(grep -h '^task dir: ' "${WORK}"/work/*/*/.command.err 2>/dev/null | sed 's/^task dir: //')
    assert_eq "SIMULATE_AND_BATCH and SCORE logged their task dir" "2" \
        "$(printf '%s\n' "${task_dirs}" | grep -c . || true)"
    assert_eq "and ran in local scratch, not the work dir" "0" \
        "$(printf '%s\n' "${task_dirs}" | grep -c "^${WORK}/work" || true)"

    # The monolithic reference: same loci, same grid, UniMap's own scripts.
    ref="${WORK}/ref"
    mkdir -p "${ref}"
    cp "${ROOT}/test/data/genome.fa" "${ref}/"
    samtools faidx "${ref}/genome.fa"
    U=(julia --startup-file=no -t 2 --project="${UNIMAP}" -m UniMap)
    (
        set -e
        cd "${ref}"
        "${U[@]}" simulate -f genome.fa -l "${WORK}/nf/loci/loci.tsv" -o sim \
            -L 30,50 -F 50,80,120,160,200 --step 5 --orientation plus
        "${U[@]}" align-script -s sim -x "${WORK}/index/genome" -t 2
        bash sim/align.sh
        "${U[@]}" score -s sim --cutoffs 0,2,5,10,20,30
    ) > "${WORK}/ref.log" 2>&1
    if [[ $? -eq 0 ]]; then ok "monolithic UniMap reference completes"; else bad "monolithic UniMap reference completes"; tail -20 "${WORK}/ref.log"; fi

    nf_scores="${WORK}/nf_scores"
    mkdir -p "${nf_scores}"
    tar -xf "${WORK}/nf/scores/mappability.zarr.tar" -C "${nf_scores}"
    cp "${WORK}/nf/scores/"*.parquet "${nf_scores}/"

    tally "$(julia --startup-file=no --project="${UNIMAP}" "${HERE}/check_scores.jl" equal "${ref}/sim/scores" "${nf_scores}" 2>&1)"
    echo
    echo "end to end -- the fixture's planted features"
    tally "$(julia --startup-file=no --project="${UNIMAP}" "${HERE}/check_scores.jl" expect "${nf_scores}" 2>&1)"

    if [[ ${DOCKER} -eq 1 ]]; then
        echo
        echo "end to end -- containers (test,docker; builds its own index)"
        # The images are amd64. On Apple Silicon they run under emulation, where
        # multithreaded Julia segfaults; one thread is slower but sound.
        JULIA_THREADS_ARG=""
        [[ "$(uname -m)" == arm64 ]] && JULIA_THREADS_ARG="--julia_threads 1"
        out=$(cd "${WORK}" && nextflow run "${ROOT}" -profile test,docker -ansi-log false \
                  --outdir "${WORK}/nfd" -w "${WORK}/workd" ${JULIA_THREADS_ARG} 2>&1)
        if [[ "${out}" == *"Done in"* ]]; then
            ok "pipeline completes in containers"
            d="${WORK}/nfd_scores"; mkdir -p "${d}"
            tar -xf "${WORK}/nfd/scores/mappability.zarr.tar" -C "${d}"
            cp "${WORK}/nfd/scores/"*.parquet "${d}/"
            tally "$(julia --startup-file=no --project="${UNIMAP}" "${HERE}/check_scores.jl" expect "${d}" 2>&1)"
        else
            bad "pipeline completes in containers"; printf '%s\n' "${out}" | tail -30
        fi
    fi
fi

# ---------------------------------------------------------------------------
echo
echo "${PASS} passed, ${FAIL} failed"
if [[ ${FAIL} -gt 0 ]]; then
    printf '  %s\n' "${FAILURES[@]}"
    exit 1
fi
