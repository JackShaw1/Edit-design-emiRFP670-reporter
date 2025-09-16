#!/bin/bash
#SBATCH --job-name=test_PRIDICT2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=24:00:00

set -euo pipefail

# ==== CONFIG ====
WORKDIR="/home/js5314/2025-09-19_edit_optimizations/PRIDICT2"
CSV_FILE="/home/js5314/2025-09-19_edit_optimizations/possible_coding_substitutions.csv"     # <-- path to your CSV with edits in col 1
# Put the original sequence here, but replace GAAT with REPLACE_ME so we can swap it per row
SEQ_TEMPLATE='CGTGTACGGTGGGAGGTCTATATAAGCAGAGCTCGTTTAGTGAACCGTCAGATCGTGTCGTGAGCTAGCTGTACACCTGCAGGAATCCCTTCTGCAGCACCTGG(+REPLACE_ME)GCCAGCTTCATCATCAGGTATTTCCGCCGAAGGCTCCGTGGCAAGACAGCCTGATCTCCTGACCTGTGAGCACGAAGAGATTCATCTGGCTGGAAGCATC'
BASE_SEQNAME="seq"

cd "$WORKDIR"

# If your CSV may have Windows newlines, normalize them on the fly.
# Also, if there’s a header, we’ll skip it by filtering out non-ACGT rows.
i=0
while IFS=, read -r EDIT _rest; do
  # skip blank lines
  [[ -z "${EDIT// }" ]] && continue
  # skip header-ish rows (anything not just A/C/G/T)
  if [[ ! "$EDIT" =~ ^[ACGTacgt]+$ ]]; then
    continue
  fi

  ((i+=1))
  # uppercase + sanitize for filenames
  EDIT_UPPER="${EDIT^^}"
  SAFE_EDIT="${EDIT_UPPER//[^A-Z0-9]/_}"

  SEQ_NAME="${BASE_SEQNAME}_${i}_${SAFE_EDIT}"
  SEQ="${SEQ_TEMPLATE//REPLACE_ME/${EDIT_UPPER}}"

  echo "[$(date +'%F %T')] Running PRIDICT2 for edit ${EDIT_UPPER} as ${SEQ_NAME} …"

  python pridict2_pegRNA_design.py single --sequence-name "${SEQ_NAME}" --sequence "${SEQ}"

  PRED_FILE="predictions/${SEQ_NAME}_pegRNA_Pridict_full.csv"

  if [[ -f "$PRED_FILE" ]]; then
    python Match_scores_coding_insertions.py "$PRED_FILE" "${EDIT_UPPER}"
    # remove the large per-edit prediction CSV after scoring (like your original script)
    # rm -f "$PRED_FILE"
  else
    echo "WARNING: Expected predictions file not found for ${SEQ_NAME}: ${PRED_FILE}" >&2
  fi

done < <(tr -d '\r' < "$CSV_FILE")

echo "[$(date +'%F %T')] All edits processed. Results are in the 'predictions' (or your scoring outputs) folder."

