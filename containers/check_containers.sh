#!/bin/bash
#SBATCH -J check-containers
#SBATCH -p cpu
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -o check_containers-%j.out
#SBATCH -e check_containers-%j.err
#
# Pull both containers and assert the commands this pipeline runs in them.
# Both a Slurm batch script and an ordinary shell script:
#
#   sbatch containers/check_containers.sh      # Unity
#   ./containers/check_containers.sh           # anywhere with docker
#
#   --runtime docker|apptainer|singularity     default: auto-detect
#   CONTAINER_UNIMAP=... CONTAINER_BOWTIE2=... override either image
#
# Keep the defaults below in step with params.container_* in nextflow.config.

set -uo pipefail

CONTAINER_UNIMAP="${CONTAINER_UNIMAP:-ghcr.io/tdw-89/unimap:0.1.0}"
CONTAINER_BOWTIE2="${CONTAINER_BOWTIE2:-quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0}"

RUNTIME=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --runtime) RUNTIME="$2"; shift 2 ;;
        -h|--help) sed -n '/^# Pull both/,/^# Keep the defaults/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *)         echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done
if [[ -z "${RUNTIME}" ]]; then
    if command -v apptainer > /dev/null 2>&1;     then RUNTIME=apptainer
    elif command -v singularity > /dev/null 2>&1; then RUNTIME=singularity
    elif command -v docker > /dev/null 2>&1;      then RUNTIME=docker
    else echo "error: found neither apptainer, singularity nor docker" >&2; exit 2
    fi
fi

# Run as the calling user, as Nextflow does -- the UniMap image must work for a
# user who cannot write to /opt.
run() {
    local image="$1"; shift
    case "${RUNTIME}" in
        docker) docker run --rm --platform linux/amd64 -u "$(id -u):$(id -g)" "${image}" bash -c "$*" ;;
        *)      "${RUNTIME}" exec "docker://${image}" bash -c "$*" ;;
    esac
}

PASS=0; FAIL=0
check() {   # check <label> <image> <command>
    local out
    if out=$(run "$2" "$3" 2>&1); then
        printf '  ok   %s\n' "$1"; PASS=$((PASS + 1))
    else
        printf '  FAIL %s\n' "$1"; printf '%s\n' "${out}" | tail -5 | sed 's/^/         /'; FAIL=$((FAIL + 1))
    fi
}

echo "runtime: ${RUNTIME}"
echo "${CONTAINER_UNIMAP}"
check "UniMap loads"                           "${CONTAINER_UNIMAP}" \
      'cd /tmp && julia --project="$UNIMAP_PROJECT" -m UniMap --version'
if [[ "${RUNTIME}" == docker ]]; then
    # Only meaningful with a fresh /tmp; Apptainer shares the host's.
    check "...from the shipped caches, compiling nothing" "${CONTAINER_UNIMAP}" \
          'cd /tmp && julia --project="$UNIMAP_PROJECT" -m UniMap --version > /dev/null && [ ! -d /tmp/julia_depot/compiled ]'
fi
check "UNIMAP_PROJECT is set and is a project" "${CONTAINER_UNIMAP}" '[ -f "$UNIMAP_PROJECT/Project.toml" ]'
check "unimap subcommands present"             "${CONTAINER_UNIMAP}" \
      'cd /tmp && julia --project="$UNIMAP_PROJECT" -m UniMap score --help | grep -q -- --fasta'
check "samtools cat / faidx"                   "${CONTAINER_UNIMAP}" 'samtools cat 2>&1 | grep -q Usage; samtools faidx 2>&1 | grep -q Usage'
check "awk, gzip, tar, xargs, cmp, git"        "${CONTAINER_UNIMAP}" 'for t in awk gzip tar xargs cmp git sort; do command -v $t > /dev/null || exit 1; done'
check "env -S shebang for check_score.jl"      "${CONTAINER_UNIMAP}" 'env -S true'
check "ps for Nextflow task metrics"           "${CONTAINER_UNIMAP}" 'command -v ps'

echo "${CONTAINER_BOWTIE2}"
check "bowtie2 and bowtie2-build"              "${CONTAINER_BOWTIE2}" 'bowtie2 --version | head -1 && bowtie2-build --version | head -1'
check "samtools view / quickcheck"             "${CONTAINER_BOWTIE2}" 'samtools --version | head -1 && samtools quickcheck 2>&1 | grep -qi usage'
check "bash with pipefail, awk"                "${CONTAINER_BOWTIE2}" 'set -o pipefail && command -v awk'

echo
echo "${PASS} passed, ${FAIL} failed"
[[ ${FAIL} -eq 0 ]]
