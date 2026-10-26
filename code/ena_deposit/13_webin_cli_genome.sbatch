#!/bin/bash
#SBATCH --job-name=ena_asm_webin
#SBATCH --output=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/assemblies/logs/webin_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/assemblies/logs/webin_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#
# Webin-CLI genome context over the assembly manifests (WORK/assemblies/manifests/*.manifest from
# 10_build_assembly_table.py), one array task per manifest. Same explicit-MODE rule as 05_webin_cli.sbatch.
#
#   MODE=validate  full local validation (sequences + metadata); nothing is sent to ENA.
#   MODE=submit    IRREVERSIBLE: validates, UPLOADS the .fa.gz over FTP, registers the analysis; ERZ returned.
#                  Keep the array width <= 12: ENA refuses FTP logins when the account holds more sessions.
#
#   N=$(ls /ibex/scratch/projects/c2014/alelopezv/ena_wgs/assemblies/manifests/*.manifest | wc -l)
#   MODE=validate sbatch --array=1-${N}%16 code/ena_deposit/13_webin_cli_genome.sbatch
#   MODE=submit   sbatch --array=1-${N}%12 code/ena_deposit/13_webin_cli_genome.sbatch     # after her go
# MANIFEST_DIR overrides the manifest folder (bins and MAGs reuse this script). Resume: aliases already OK in the
# status log exit immediately. Credential: ~/.webin_pass via -passwordFile (never argv).

set -u
WORK=/ibex/scratch/projects/c2014/alelopezv/ena_wgs
MANIFEST_DIR=${MANIFEST_DIR:-${WORK}/assemblies/manifests}
LOGDIR=${LOGDIR:-${WORK}/assemblies/logs}
WEBIN_CLI=${WEBIN_CLI:-${WORK}/webin-cli-10.0.0.jar}
WEBIN_USER=Webin-65194
PASS_FILE=${HOME}/.webin_pass
MODE=${MODE:-}

case "${MODE}" in
    validate|submit) ;;
    *) echo "ERROR: set MODE=validate (dry run) or MODE=submit (irreversible). No default."; exit 1 ;;
esac
[ -s "${PASS_FILE}" ] || { echo "ERROR: ${PASS_FILE} missing/empty"; exit 1; }
[ -f "${WEBIN_CLI}" ] || { echo "ERROR: ${WEBIN_CLI} missing"; exit 1; }
module load java/19.0.1 || { echo "ERROR: java/19.0.1 module not available (Webin-CLI needs Java 17+)"; exit 1; }

STATUS=${LOGDIR}/genome_${MODE}_status.tsv     # alias  status  erz  timestamp
REPORTS=${LOGDIR}/webin_reports_${MODE}
mkdir -p "${REPORTS}"
[ -f "${STATUS}" ] || printf 'alias\tstatus\terz\ttimestamp\n' > "${STATUS}"

mapfile -t manifests < <(ls "${MANIFEST_DIR}"/*.manifest | sort)
manifest=${manifests[$(( SLURM_ARRAY_TASK_ID - 1 ))]:-}
[ -n "${manifest}" ] || { echo "ERROR: no manifest for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }
alias=$(basename "${manifest}" .manifest)

if grep -q "^${alias}	OK	" "${STATUS}"; then
    echo "[skip] ${alias} already ${MODE} OK"; exit 0
fi

echo "[${MODE}] ${alias} $(date -Is)"
output=$(java -jar "${WEBIN_CLI}" -context genome -manifest "${manifest}" \
            -userName "${WEBIN_USER}" -passwordFile "${PASS_FILE}" \
            -${MODE} -outputDir "${REPORTS}" 2>&1)
rc=$?
echo "${output}" | tail -6
if [ ${rc} -eq 0 ]; then
    erz=$(echo "${output}" | grep -o 'ERZ[0-9]*' | tail -1)
    printf '%s\tOK\t%s\t%s\n' "${alias}" "${erz}" "$(date -Is)" >> "${STATUS}"
else
    printf '%s\tFAILED\t\t%s\n' "${alias}" "$(date -Is)" >> "${STATUS}"
    echo "  FAILED - see ${REPORTS}/genome/${alias}/"
fi
exit ${rc}
