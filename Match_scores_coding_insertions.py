import sys
import os
import csv

def reverse_and_switch(seq):
    comp = {'A':'T','T':'A','G':'C','C':'G','a':'t','t':'a','g':'c','c':'g'}
    return ''.join(comp[b] for b in seq[::-1])

def best_substring_similarity(target, text):
    """
    Return the best fraction of matching characters between target
    and any same-length window in text (no indels).
    """
    target = target.upper()
    text = text.upper()
    n, m = len(target), len(text)
    if n == 0 or m == 0 or m < n:
        return 0.0
    best = 0
    for i in range(m - n + 1):
        window = text[i:i+n]
        matches = sum(1 for a,b in zip(target, window) if a == b)
        frac = matches / n
        if frac > best:
            best = frac
            if best == 1.0:
                break
    return best

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <predictions_csv> <INSERT_SEQ>")
        sys.exit(1)

    csv_path = sys.argv[1]
    insert_seq = sys.argv[2].strip().upper()

    coding_edit_spacer = 'GGAATCCCTTCTGCAGCACC'  # 20 nt
    # RTT
    new_seq_rc = reverse_and_switch(insert_seq).upper()
    coding_edit_rtt   = f'GATGATGAAGCTGGC{new_seq_rc}CCAGGT'
    coding_edit_pbs   = 'GCTGCAGAAGGGATT'

    # weights for finding most similar row
    w_rtt, w_spacer, w_pbs = 0.5, 0.25, 0.25

    best_row = None
    best_score = -1.0

    with open(csv_path, 'r', newline='') as f:
        reader = csv.reader(f)
        for row in reader:

            if len(row) < 24:
                continue

            pot_spacer = (row[20] or "").strip()
            pot_pbs    = (row[21] or "").strip()
            pot_rtt    = (row[23] or "").strip()

            # per-feature similarity (0..1), using best substring match
            s_rtt    = best_substring_similarity(coding_edit_rtt, pot_rtt)
            s_spacer = best_substring_similarity(coding_edit_spacer, pot_spacer)
            s_pbs    = best_substring_similarity(coding_edit_pbs, pot_pbs)

            # combine
            score = w_rtt * s_rtt + w_spacer * s_spacer + w_pbs * s_pbs

            if score > best_score:
                best_score = score
                best_row = row

    if best_row is None:
        print("No usable rows found.")
        sys.exit(2)


    out_path = 'coding_insert_scores_sep_16.csv'
    with open(out_path, 'a', newline='') as out:
        writer = csv.writer(out)
        writer.writerow([insert_seq, best_row[1], f"{best_score:.4f}"])

    print(f"Best row chosen with combined similarity {best_score:.4f}. Wrote to {out_path}.")

if __name__ == "__main__":
    main()
