#!/bin/bash
#SBATCH --job-name=ena_wgs_upload
#SBATCH --output=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/logs/upload_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/logs/upload_%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#
# md5 + FTP upload of the WGS read files to the ENA Webin upload area, one Slurm array task per run
# (R1 + R2). Adapted from eqmgc 00_upload_ftp.sh (2026-04-30) with the IFS and log-overwrite fixes kept.
#
# Reads the symlink farm built by 00_build_runs_table.py --symlinks; uploads under the clean public
# name (basename of the symlink), so nothing on disk is renamed or copied.
#
# Usage (array size = number of data rows in wgs_runs.tsv; %8 = concurrent uploads):
#   N=$(($(wc -l < code/ena_deposit/wgs_runs.tsv) - 1))
#   sbatch --array=1-${N}%8 code/ena_deposit/01_upload_ftp.sbatch
# Re-running skips files already logged OK (resume). DRYRUN=1 computes md5 only, uploads nothing.
#
# Credential: ~/.webin_pass (mode 600), passed to curl as a config file - never on the command line.

set -u
WORK=/ibex/scratch/projects/c2014/alelopezv/ena_wgs
TSV=/ibex/scratch/projects/c2014/alelopezv/repos/empty-quarter-mags/code/ena_deposit/wgs_runs.tsv
WEBIN_HOST=webin2.ebi.ac.uk
PASS_FILE=${HOME}/.webin_pass
MD5LOG=${WORK}/logs/md5sums.tsv          # upload_name  md5  bytes
UPLOG=${WORK}/logs/upload_status.tsv     # upload_name  status  remote_bytes  timestamp
DRYRUN=${DRYRUN:-0}

[ -s "${PASS_FILE}" ] || { echo "ERROR: ${PASS_FILE} missing/empty"; exit 1; }
[ -f "${MD5LOG}" ] || printf 'upload_name\tmd5\tbytes\n' > "${MD5LOG}"
[ -f "${UPLOG}" ]  || printf 'upload_name\tstatus\tremote_bytes\ttimestamp\n' > "${UPLOG}"

# row = SLURM_ARRAY_TASK_ID (1-based data row); columns resolved by header name
row=$(( SLURM_ARRAY_TASK_ID + 1 ))
line=$(sed -n "${row}p" "${TSV}")
[ -n "${line}" ] || { echo "ERROR: no row ${row} in ${TSV}"; exit 1; }
col() { head -1 "${TSV}" | tr '\t' '\n' | grep -n -x "$1" | cut -d: -f1; }
c_alias=$(col run_alias); c_r1=$(col upload_r1); c_r2=$(col upload_r2)
[ -n "${c_alias}" ] && [ -n "${c_r1}" ] && [ -n "${c_r2}" ] || { echo "ERROR: columns run_alias/upload_r1/upload_r2 not found in ${TSV}"; exit 1; }
run_alias=$(printf '%s' "${line}" | cut -f"${c_alias}")
files=( "$(printf '%s' "${line}" | cut -f"${c_r1}")" "$(printf '%s' "${line}" | cut -f"${c_r2}")" )

curlconf=$(mktemp)
trap 'rm -f "${curlconf}"' EXIT
printf 'user = "Webin-65194:%s"\n' "$(cat "${PASS_FILE}")" > "${curlconf}"

rc_all=0
for name in "${files[@]}"; do
    path="${WORK}/fastq/${name}"
    [ -e "${path}" ] || { echo "ERROR: ${path} missing"; rc_all=1; continue; }
    bytes=$(stat -L -c %s "${path}")

    if ! grep -q "^${name}	" "${MD5LOG}"; then
        md5=$(md5sum "${path}" | cut -d' ' -f1)
        printf '%s\t%s\t%s\n' "${name}" "${md5}" "${bytes}" >> "${MD5LOG}"
        echo "[md5] ${name} ${md5} ${bytes}"
    else
        echo "[md5] ${name} already logged"
    fi

    [ "${DRYRUN}" = "1" ] && { echo "[dryrun] would upload ${name}"; continue; }
    if grep -q "^${name}	OK	" "${UPLOG}"; then echo "[skip] ${name} already uploaded"; continue; fi

    echo "[upload] ${run_alias} ${name} ($(( bytes / 1000000 )) MB) ..."
    if curl --silent --show-error -K "${curlconf}" --ftp-pasv -T "${path}" "ftp://${WEBIN_HOST}/${name}"; then
        remote=$(curl --silent -K "${curlconf}" --ftp-pasv -I "ftp://${WEBIN_HOST}/${name}" | awk '/Content-Length/ {print $2}' | tr -d '\r')
        if [ "${remote}" = "${bytes}" ]; then
            printf '%s\tOK\t%s\t%s\n' "${name}" "${remote}" "$(date -Is)" >> "${UPLOG}"
            echo "  -> OK (remote ${remote} bytes)"
        else
            printf '%s\tSIZE_MISMATCH\t%s\t%s\n' "${name}" "${remote:-?}" "$(date -Is)" >> "${UPLOG}"
            echo "  -> SIZE MISMATCH local ${bytes} remote ${remote:-?}"; rc_all=1
        fi
    else
        printf '%s\tFAILED\t\t%s\n' "${name}" "$(date -Is)" >> "${UPLOG}"
        echo "  -> FAILED"; rc_all=1
    fi
done
exit ${rc_all}
