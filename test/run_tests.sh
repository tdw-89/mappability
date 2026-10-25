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

"${BIN}/plan_batches.sh" regions.tsv tiling.tsv 2 plan.tsv batches.tsv 2> /dev/null
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

"${BIN}/plan_batches.sh" regions.tsv tiling.tsv 40 plan40.tsv batches40.tsv 2> /dev/null
assert_eq "n_batches capped at the region count" "5" "$(tail -n +2 plan40.tsv | cut -f2 | sort -u | wc -l | tr -d ' ')"

grep -v '^r3' tiling.tsv > tiling_missing.tsv
"${BIN}/plan_batches.sh" regions.tsv tiling_missing.tsv 2 p.tsv b.tsv 2> err.txt
rc=$?
if [[ ${rc} -ne 0 ]]; then ok "region absent from tiling.tsv is an error"; else bad "region absent from tiling.tsv is an error"; fi
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
out=$(nf_launch --fragment_lengths '50,abc' -stub-run)
assert_contains "malformed --fragment_lengths is refused" "comma-separated list" "${out}"
out=$(nf_launch --bowtie2_index /nonexistent/idx -stub-run)
assert_contains "incomplete --bowtie2_index is refused" "no complete bowtie2 index" "${out}"
out=$(cd "${WORK}" && nextflow run "${ROOT}" -profile test -ansi-log false -stub-run 2>&1)
assert_contains "no container and no --unimap_project is refused" "--unimap_project" "${out}"
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

    out=$(nf_launch --bowtie2_index "${WORK}/index/genome" --outdir "${WORK}/nf" -w "${WORK}/work")
    if [[ "${out}" == *"Done in"* ]]; then
        ok "pipeline completes"
    else
        bad "pipeline completes"; printf '%s\n' "${out}" | tail -30
    fi

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
