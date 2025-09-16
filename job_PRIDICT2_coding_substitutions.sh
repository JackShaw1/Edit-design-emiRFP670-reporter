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

# IMPORTANT: Keep a single explicit marker for the canonical 4-nt block (CTGG).
# We'll adapt the replaced span to the edit length, but we ALWAYS output one "(REF/ALT)" pair.
SEQ_TEMPLATE='CGTGTACGGTGGGAGGTCTATATAAGCAGAGCTCGTTTAGTGAACCGTCAGATCGTGTCGTGAGCTAGCTGTACACCTGCAGGAATCCCTTCTGCAGCAC(CTGG/REPLACE_ME)GCCAGCTTCATCATCAGGTATTTCCGCCGAAGGCTCCGTGGCAAGACAGCCTGATCTCCTGACCTGTGAGCACGAAGAGATTCATCTGGCTGGAAGCATC'
BASE_SEQNAME="seq"

cd "$WORKDIR"

MARK="(CTGG/REPLACE_ME)"
ORIG_BLOCK="CTGG"
ORIG_LEN=${#ORIG_BLOCK}

# Split template around the MARK once
if [[ "$SEQ_TEMPLATE" != *"$MARK"* ]]; then
  echo "ERROR: SEQ_TEMPLATE must contain the marker $MARK exactly once." >&2
  exit 1
fi
LEFT="${SEQ_TEMPLATE%%$MARK*}"
RIGHT="${SEQ_TEMPLATE#*"$MARK"}"

i=0
while IFS=, read -r EDIT _rest; do
  # skip blank lines
  [[ -z "${EDIT// }" ]] && continue
  # skip header-ish rows (anything not just A/C/G/T)
  if [[ ! "$EDIT" =~ ^[ACGTacgt]+$ ]]; then
    continue
  fi

  ((i+=1))
  EDIT_UPPER="${EDIT^^}"
  SAFE_EDIT="${EDIT_UPPER//[^A-Z0-9]/_}"
  SEQ_NAME="${BASE_SEQNAME}_${i}_${SAFE_EDIT}"

  EDIT_LEN=${#EDIT_UPPER}

  # - If EDIT_LEN < 4:
  #     keep the leftmost (4 - EDIT_LEN) bases of CTGG outside the brackets,
  #     and replace the remainder of CTGG INSIDE the brackets as (REF/ALT).
  #     Example: EDIT_LEN=2 -> CT(GG/NN)
  # - If EDIT_LEN = 4:
  #     replace CTGG entirely: (CTGG/NNNN)
  # - If EDIT_LEN > 4:
  #     expand the REF span to the LEFT by (EDIT_LEN-4) bases upstream of CTGG,
  #     so REF = upstream_span + CTGG, e.g. EDIT_LEN=6 -> (ACCTGG/NNNNNN)
  if (( EDIT_LEN < ORIG_LEN )); then
    KEEP_LEFT_COUNT=$(( ORIG_LEN - EDIT_LEN ))
    KEPT_OUTSIDE="${ORIG_BLOCK:0:KEEP_LEFT_COUNT}"
    REF_SPAN="${ORIG_BLOCK:KEEP_LEFT_COUNT}"
    # Final sequence: LEFT + kept outside + (REF/ALT) + RIGHT
    SEQ="${LEFT}${KEPT_OUTSIDE}(${REF_SPAN}/${EDIT_UPPER})${RIGHT}"

  elif (( EDIT_LEN == ORIG_LEN )); then
    REF_SPAN="${ORIG_BLOCK}"
    SEQ="${LEFT}(${REF_SPAN}/${EDIT_UPPER})${RIGHT}"

  else
    EXTRA_LEFT=$(( EDIT_LEN - ORIG_LEN ))
    LEFT_LEN=${#LEFT}
    if (( EXTRA_LEFT > LEFT_LEN )); then
      echo "WARNING: edit '${EDIT_UPPER}' needs ${EXTRA_LEFT} upstream bases but LEFT has only ${LEFT_LEN}. Skipping." >&2
      continue
    fi
    # Upstream bases to include in REF (and remove from LEFT outside):
    UPSTREAM_SPAN="${LEFT: -$EXTRA_LEFT}"
    NEW_LEFT="${LEFT:0:$(( LEFT_LEN - EXTRA_LEFT ))}"
    REF_SPAN="${UPSTREAM_SPAN}${ORIG_BLOCK}"
    SEQ="${NEW_LEFT}(${REF_SPAN}/${EDIT_UPPER})${RIGHT}"
  fi

  echo "[$(date +'%F %T')] Running PRIDICT2 for edit ${EDIT_UPPER} as ${SEQ_NAME} …"

  # NOTE: SEQ now contains exactly one '(REF/ALT)' pair as required by PRIDICT2.
  python pridict2_pegRNA_design.py single --sequence-name "${SEQ_NAME}" --sequence "${SEQ}"

  PRED_FILE="predictions/${SEQ_NAME}_pegRNA_Pridict_full.csv"

  if [[ -f "$PRED_FILE" ]]; then
    python Match_scores_coding_subs.py "$PRED_FILE" "${EDIT_UPPER}"
    # rm -f "$PRED_FILE"  # optional cleanup
  else
    echo "WARNING: Expected predictions file not found for ${SEQ_NAME}: ${PRED_FILE}" >&2
  fi

done < <(tr -d '\r' < "$CSV_FILE")

echo "[$(date +'%F %T')] All edits processed. Results are in the 'predictions' (and your scoring outputs) folder."
