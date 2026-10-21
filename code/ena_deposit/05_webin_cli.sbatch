#!/bin/bash
#SBATCH --job-name=ena_wgs_webin
#SBATCH --output=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/logs/webin_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/logs/webin_%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#
# Webin-CLI over the reads manifests, one array task per manifest (validation reads the whole FASTQ
# pair, so 554 runs must not be serialised in one job). Merged from eqmgc 02_submit_webin.sh +
# 03_sbatch_submit.sh; the explicit-MODE rule is kept: no default, no accidental submit.
#
#   MODE=validate  local check of manifest + files (nothing is sent to ENA; the files need not be
#                  uploaded yet) - this is the dry run.
#   MODE=submit    IRREVERSIBLE: registers experiment + run on ENA. ⚠ Webin-CLI UPLOADS THE FILES ITSELF into
#                  webin-cli/reads/<run>/ and ignores files pre-uploaded to the root of the upload area
#                  (learned 2026-09-06 after 46 runs): for pre-uploaded reads use 08_submit_runs_xml.py instead.
#   Never combine -submit with -test: Webin-CLI would push the 6 TB to the test area. Not offered here.
#
# Usage:
#   N=$(ls /ibex/scratch/projects/c2014/alelopezv/ena_wgs/manifests/*.manifest | wc -l)
#   MODE=validate sbatch --array=1-${N}%16 code/ena_deposit/05_webin_cli.sbatch
#   MODE=submit   sbatch --array=1-${N}%4  code/ena_deposit/05_webin_cli.sbatch     # after her go
# Resume: tasks whose alias is already OK in the status log exit immediately.
# Credential: ~/.webin_pass, read into a variable and passed with -passwordFile (never argv).

set -u
WORK=/ibex/scratch/projects/c2014/alelopezv/ena_wgs
MANIFEST_DIR=${MANIFEST_DIR:-${WORK}/manifests}
WEBIN_CLI=${WORK}/webin-cli.jar
WEBIN_USER=Webin-65194
PASS_FILE=${HOME}/.webin_pass
MODE=${MODE:-}

case "${MODE}" in
    validate|submit) ;;
    *) echo "ERROR: set MODE=validate (dry run) or MODE=submit (irreversible). No default."; exit 1 ;;
esac
[ -s "${PASS_FILE}" ] || { echo "ERROR: ${PASS_FILE} missing/empty"; exit 1; }
[ -f "${WEBIN_CLI}" ] || { echo "ERROR: ${WEBIN_CLI} missing (cp from eqmgc/webin-cli.jar, v9.0.3)"; exit 1; }
module load java/19.0.1 || { echo "ERROR: java/19.0.1 module not available (Webin-CLI needs Java 17+)"; exit 1; }

STATUS=${WORK}/logs/webin_${MODE}_status.tsv     # alias  status  ERX  ERR  timestamp
REPORTS=${WORK}/logs/webin_reports_${MODE}
mkdir -p "${REPORTS}"
[ -f "${STATUS}" ] || printf 'alias\tstatus\terx\terr\ttimestamp\n' > "${STATUS}"

mapfile -t manifests < <(ls "${MANIFEST_DIR}"/*.manifest | sort)
manifest=${manifests[$(( SLURM_ARRAY_TASK_ID - 1 ))]:-}
[ -n "${manifest}" ] || { echo "ERROR: no manifest for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }
alias=$(basename "${manifest}" .manifest)

if [ "${MODE}" = "submit" ] && grep -q "^${alias}	OK	" "${STATUS}"; then
    echo "[skip] ${alias} already submitted OK"; exit 0
fi

echo "[${MODE}] ${alias}"
output=$(java -jar "${WEBIN_CLI}" -context reads -manifest "${manifest}" \
            -userName "${WEBIN_USER}" -passwordFile "${PASS_FILE}" \
            -${MODE} -outputDir "${REPORTS}" 2>&1)
rc=$?
echo "${output}" | tail -5
if [ ${rc} -eq 0 ]; then
    erx=$(echo "${output}" | grep -o 'ERX[0-9]*' | tail -1)
    err=$(echo "${output}" | grep -o 'ERR[0-9]*' | tail -1)
    printf '%s\tOK\t%s\t%s\t%s\n' "${alias}" "${erx}" "${err}" "$(date -Is)" >> "${STATUS}"
else
    printf '%s\tFAILED\t\t\t%s\n' "${alias}" "$(date -Is)" >> "${STATUS}"
    echo "  FAILED - see ${REPORTS}/${alias}/"
fi
exit ${rc}
