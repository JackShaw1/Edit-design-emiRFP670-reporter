#!/bin/bash
#SBATCH --job-name=test_PRIDICT2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

# ==== CONFIG ====
WORKDIR="/home/js5314/2025-09-19_edit_optimizations/PRIDICT2"
CSV_FILE="/home/js5314/2025-09-19_edit_optimizations/possible_template_substitutions.csv"

# Anchor the canonical 4-nt block CTGG where the edit lands.
# We'll adapt the replaced span based on the edit length and ALWAYS emit exactly one "(REF/ALT)" block.
SEQ_TEMPLATE='GGCTGGATGCTTCCAGCCAGATGAATCTCTTCGTGCTCACAGGTCAGGAGATCAGGCTGTCTTGCCACGGAGCCTTCGGCGGAAATACCTGATGATGAAG(CTGG/REPLACE_ME)CCCAGGTGCTGCAGAAGGGATTCCTGCAGGTGTACAGCTAGCTCACGACACGATCTGACGGTTCACTAAACGAGCTCTGCTTATATAGACCTCCCACCGT'
BASE_SEQNAME="seq"

cd "$WORKDIR"

MARK="(CTGG/REPLACE_ME)"
ORIG_BLOCK="CTGG"
ORIG_LEN=${#ORIG_BLOCK}

# Split template once around the marker
if [[ "$SEQ_TEMPLATE" != *"$MARK"* ]]; then
  echo "ERROR: SEQ_TEMPLATE must contain the marker $MARK exactly once." >&2
  exit 1
fi
LEFT="${SEQ_TEMPLATE%%$MARK*}"
RIGHT="${SEQ_TEMPLATE#*"$MARK"}"

# -------------------- helpers --------------------
substr() { local s="$1" pos="$2" len="$3"; echo -n "${s:pos:len}"; }

common_prefix_len() {
  local r="$1" a="$2"
  local lr=${#r} la=${#a}
  local m=$(( lr < la ? lr : la ))
  local i=0
  while (( i < m )) && [[ "${r:i:1}" == "${a:i:1}" ]]; do ((i++)); done
  echo -n "$i"
}

common_suffix_len() {
  local r="$1" a="$2"
  local lr=${#r} la=${#a}
  local m=$(( lr < la ? lr : la ))
  local i=0
  while (( i < m )) && [[ "${r:lr-1-i:1}" == "${a:la-1-i:1}" ]]; do ((i++)); done
  echo -n "$i"
}
# -------------------------------------------------

i=0
while IFS=, read -r EDIT _rest; do
  # skip blank lines / headers
  [[ -z "${EDIT// }" ]] && continue
  if [[ ! "$EDIT" =~ ^[ACGTacgt]+$ ]]; then
    continue
  fi

  ((i+=1))
  EDIT_UPPER="${EDIT^^}"
  SAFE_EDIT="${EDIT_UPPER//[^A-Z0-9]/_}"
  SEQ_NAME="${BASE_SEQNAME}_${i}_${SAFE_EDIT}"

  EDIT_LEN=${#EDIT_UPPER}

  # Choose REF_SPAN based on edit length (relative to CTGG)
  REF_SPAN=""
  NEW_LEFT="$LEFT"
  NEW_RIGHT="$RIGHT"

  if (( EDIT_LEN < ORIG_LEN )); then
    # keep leftmost (4 - L) of CTGG outside; replace the remainder inside ()
    KEEP_LEFT_COUNT=$(( ORIG_LEN - EDIT_LEN ))
    KEPT_OUTSIDE="${ORIG_BLOCK:0:KEEP_LEFT_COUNT}"
    REF_SPAN="${ORIG_BLOCK:KEEP_LEFT_COUNT}"
    NEW_LEFT="${LEFT}${KEPT_OUTSIDE}"

  elif (( EDIT_LEN == ORIG_LEN )); then
    # replace CTGG entirely
    REF_SPAN="${ORIG_BLOCK}"

  else
    # EDIT_LEN > 4 → expand replacement to the LEFT by (EDIT_LEN - 4) upstream bases
    EXTRA_LEFT=$(( EDIT_LEN - ORIG_LEN ))
    LEFT_LEN=${#LEFT}
    if (( EXTRA_LEFT > LEFT_LEN )); then
      echo "WARNING: '${EDIT_UPPER}' needs ${EXTRA_LEFT} upstream bases but LEFT has ${LEFT_LEN}. Skipping ${SEQ_NAME}." >&2
      continue
    fi
    UPSTREAM_SPAN="${LEFT: -$EXTRA_LEFT}"
    NEW_LEFT="${LEFT:0:$(( LEFT_LEN - EXTRA_LEFT ))}"
    REF_SPAN="${UPSTREAM_SPAN}${ORIG_BLOCK}"
  fi

  # Trim shared prefix/suffix so only the true diff remains in (REF/ALT)
  PREFIX=$(common_prefix_len "$REF_SPAN" "$EDIT_UPPER")
  REF_CORE="${REF_SPAN:PREFIX}"
  ALT_CORE="${EDIT_UPPER:PREFIX}"
  SUFFIX=$(common_suffix_len "$REF_CORE" "$ALT_CORE")

  if (( PREFIX > 0 )); then
    NEW_LEFT="${NEW_LEFT}$(substr "$REF_SPAN" 0 "$PREFIX")"
  fi
  if (( SUFFIX > 0 )); then
    NEW_RIGHT="$(substr "$REF_CORE" $(( ${#REF_CORE} - SUFFIX )) "$SUFFIX")${NEW_RIGHT}"
  fi

  REF_CORE="$(substr "$REF_CORE" 0 $(( ${#REF_CORE} - SUFFIX )))"
  ALT_CORE="$(substr "$ALT_CORE" 0 $(( ${#ALT_CORE} - SUFFIX )))"

  # If nothing differs after trimming, skip
  if [[ -z "$REF_CORE" && -z "$ALT_CORE" ]]; then
    echo "WARNING: Edit '${EDIT_UPPER}' identical to reference after trimming; skipping ${SEQ_NAME}." >&2
    continue
  fi

  # Build final sequence with EXACTLY ONE (REF/ALT)
  SEQ="${NEW_LEFT}(${REF_CORE}/${ALT_CORE})${NEW_RIGHT}"

  # sanity: exactly one '(' and ')'
  if [[ "$(grep -o "(" <<<"$SEQ" | wc -l)" -ne 1 || "$(grep -o ")" <<<"$SEQ" | wc -l)" -ne 1 ]]; then
    echo "ERROR: Generated sequence has != 1 bracket pair for ${SEQ_NAME}. Skipping." >&2
    continue
  fi

  echo "[$(date +'%F %T')] Running PRIDICT2 for edit ${EDIT_UPPER} as ${SEQ_NAME} …"
  python pridict2_pegRNA_design_template.py single --sequence-name "${SEQ_NAME}" --sequence "${SEQ}"

  PRED_FILE="predictions/${SEQ_NAME}_pegRNA_Pridict_full_template.csv"
  if [[ -f "$PRED_FILE" ]]; then
    python Match_scores_template_subs.py "$PRED_FILE" "${EDIT_UPPER}"
    # rm -f "$PRED_FILE"   # optional cleanup
  else
    echo "WARNING: Expected predictions file not found for ${SEQ_NAME}: ${PRED_FILE}" >&2
  fi

done < <(tr -d '\r' < "$CSV_FILE")

echo "[$(date +'%F %T')] All edits processed. Results are in the 'predictions' (or your scoring outputs) folder."
