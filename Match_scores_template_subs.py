import sys
import os
sys.path.insert(0, os.path.expanduser("PRIDICT2/python_libraries"))
import csv

def reverse_and_switch(seq):
    dict_chars = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }
    seq = seq[::-1].upper()
    new_seq = ''
    for char in seq:
        new_char = dict_chars[char]
        new_seq = new_seq + new_char
    return new_seq

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <pred_csv> <edit>")
        sys.exit(1)

    pred_csv = sys.argv[1]
    edit = sys.argv[2].upper()

    new_seq = reverse_and_switch(edit).upper()


    coding_edit_spacer = 'GAAATACCTGATGATGAAGC'
    double = False
    if len(edit) == 3:
        coding_edit_rtt    = reverse_and_switch(f'AGC{edit}CCCAGGTGC').upper()
    if len(edit) == 4:
        coding_edit_rtt    = reverse_and_switch(f'AG{edit}CCCAGGTGC').upper()
    if len(edit) == 5:
        coding_edit_rtt    = reverse_and_switch(f'A{edit}CCCAGGTGC').upper()
    if len(edit) == 6:
        coding_edit_rtt    = reverse_and_switch(f'{edit}CCCAGGTGC').upper()
    coding_edit_pbs    = 'TCATCATCAGGTATT'


    taken = []
    wrote = False
    with open(pred_csv, 'r', newline='') as file, \
        open('template_sub_scores_sep_15.csv', 'a', newline='') as outfile:

        reader = csv.reader(file)
        writer = csv.writer(outfile)

        for row in reader:

            pot_spacer = row[20].upper()
            pot_pbs    = row[21].upper()
            pot_rtt    = row[23].upper()

            if (coding_edit_rtt in pot_rtt and coding_edit_pbs in pot_pbs and coding_edit_spacer in pot_spacer and coding_edit_rtt not in taken):
                writer.writerow([edit, row[1]])
                wrote = True
                taken.append(coding_edit_rtt)



    if not wrote:
        print(f"[mapper] No matching rows found for edit {edit}. "
              f"Check column indices (20/21/23) and strings.", file=sys.stderr)




if __name__ == "__main__":

    print('HELLO WORLD')
    main()