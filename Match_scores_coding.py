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
        print(f"Usage: {sys.argv[0]} <argument>")
        sys.exit(1)
    new_seq = reverse_and_switch(sys.argv[2]).lower()
    with open(sys.argv[1], 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            identity = False
            pot_spacer = row[20]
            pot_rtt = row[23]
            pot_pbs = row[21]
            coding_edit_spacer = 'GGAATCCCTTCTGCAGCACC'
            if len(new_seq) == 6:
                coding_edit_rtt = f'GATGATGAAGCTGGC{new_seq}'
            if len(new_seq) == 5:
                coding_edit_rtt = f'GATGATGAAGCTGGC{new_seq}T'
            if len(new_seq) == 4:
                coding_edit_rtt = f'GATGATGAAGCTGGC{new_seq}GT'
            if len(new_seq) == 3:
                coding_edit_rtt = f'GATGATGAAGCTGGC{new_seq}GGT'
            if len(new_seq) == 2:
                coding_edit_rtt = f'GATGATGAAGCTGGC{new_seq}AGGT'
            coding_edit_pbs = 'GCTGCAGAAGGGATT'
            if coding_edit_rtt == pot_rtt:
                if coding_edit_pbs == pot_pbs:
                    if coding_edit_spacer == pot_spacer:
                        with open('coding_sub_scores_sep_15.csv', 'a', newline='') as outfile:
                            writer = csv.writer(outfile)
                            writer.writerow([sys.argv[2], row[1]])
                        identity = True
            if identity == True:
                print(f'identity: {identity}')



if __name__ == "__main__":
    main()