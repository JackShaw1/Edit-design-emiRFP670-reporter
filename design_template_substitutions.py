import itertools
import csv
import random
import math

genetic_code = {
    # Phenylalanine / Leucine
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    # Leucine
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    # Isoleucine / Methionine (start)
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    # Valine
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    # Serine
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "AGT":"S","AGC":"S",
    # Proline
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    # Threonine
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    # Alanine
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    # Tyrosine / Stop
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    # Histidine / Glutamine
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    # Asparagine / Lysine
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    # Aspartic acid / Glutamic acid
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    # Cysteine / Tryptophan / Stop
    "TGT":"C","TGC":"C","TGG":"W","TGA":"*",
    # Arginine
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGA":"R","AGG":"R",
    # Glycine
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

nucleotides = ['A', 'T', 'G', 'C']

# Create and combine all possible substitution edits into one long list 
all_possible_edits_3 = [''.join(p) for p in itertools.product(nucleotides, repeat=3)]
all_possible_edits_4 = [''.join(p) for p in itertools.product(nucleotides, repeat=4)]
all_possible_edits_5 = [''.join(p) for p in itertools.product(nucleotides, repeat=5)]
all_possible_edits_6 = [''.join(p) for p in itertools.product(nucleotides, repeat=6)]

all_possible_edits = all_possible_edits_3 + all_possible_edits_4 + all_possible_edits_5 + all_possible_edits_6

start_pos_dict = {
    3: 4,
    4: 3,
    5: 2,
    6: 1
}

def calc_gc_content(seq):
    total = 0
    gc = 0
    edit = seq.upper()
    for char in edit:
        total += 1
        if char in ['G', 'C']:
            gc += 1
    return gc/total

def reverse_complement(seq: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join(complement[base] for base in reversed(seq))

# For loop iterates over all possible substitutions at our locus and gathers metrics for them
counter = 0
OG_seq = 'AGCTGG'
for possible_edit in all_possible_edits:

    # Find start codon
    if possible_edit[-3:] == 'CAT':
        
        new_seq = OG_seq[:start_pos_dict[len(possible_edit)]] + possible_edit
        identity = True

        # No redundant sequences
        if len(new_seq) == 6 and new_seq[:4] == 'AGCT':
            identity = False
        if len(new_seq) == 5 and new_seq[:3] == 'GCT':
            identity = False
        if len(new_seq) == 4 and new_seq[:2] == 'CT':
            identity = False
        if len(new_seq) == 3 and new_seq[:1] == 'T':
            identity = False

        
        if identity == True:
            new_reporter_seq = None

            # Ensure no out of frame edits
            if len(possible_edit) == 3:
                new_reporter_seq = 'MASSSGISAEG'
            else:
                if genetic_code[reverse_complement(possible_edit[:3])] == '*':
                    continue
                else:
                    new_reporter_seq = 'M' + genetic_code[reverse_complement(possible_edit[:3])] + 'PASSSGISAEG'

            new_rtt = 'GCACCTGGGATGCCT' + reverse_complement(possible_edit)
            gc_content = calc_gc_content(possible_edit)
            length = len(possible_edit)
            position = start_pos_dict[len(possible_edit)]
            with open('possible_template_substitutions.csv', 'a', newline='') as outfile:
                writer = csv.writer(outfile)
                writer.writerow([possible_edit, OG_seq, new_seq, position, gc_content, length, new_rtt, new_reporter_seq])
