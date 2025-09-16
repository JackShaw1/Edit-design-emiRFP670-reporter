import itertools
import csv

def translate_from_start(seq: str) -> str:
    """
    Translate a DNA sequence into amino acids starting at the first ATG.
    Returns the amino acid sequence as a string of one-letter codes.
    Stops at the first stop codon (*).
    """
    # standard codon table
    codon_table = {
        # Phenylalanine / Leucine
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
        "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
        "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
        "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
        "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
        "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
        "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
        "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
        "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
        "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
        "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
        "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
        "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
        "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
        "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
    }

    seq = seq.upper()
    
    # find the first ATG
    start_index = seq.find("ATG")
    if start_index == -1:
        return ""   # no start codon

    protein = []
    # step through codons from the start
    for i in range(start_index, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            break
        aa = codon_table.get(codon, "X")
        if aa == "*":   # stop codon
            break
        protein.append(aa)

    return "".join(protein)

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

def generate_dna_sequences(max_len):
    bases = ['A', 'C', 'G', 'T']
    all_seqs = []
    for length in range(1, max_len + 1):
        for combo in itertools.product(bases, repeat=length):
            all_seqs.append(''.join(combo))
    return all_seqs

original_seq = 'ACctggGccagCTTCATCATCAGGTATTTCC'.upper()

max_length = 6
POSITION = 1

all_possible_seqs = generate_dna_sequences(max_length)

for seq in all_possible_seqs:
    if seq[-3:] != 'ATG' and seq[-4:-1] != 'ATG':

        identity = False

        # Check for valid insertions
        if len(seq) == 6 and seq[1:4] == 'ATG':
            identity = True
        if len(seq) == 5 and seq[0:3] == 'ATG':
            identity = True

        if identity == True:
            
            total_seq = seq + original_seq
            RFP670_seq = translate_from_start(total_seq)
            gc_content = calc_gc_content(seq)
            length = len(seq)

            if len(RFP670_seq) > 4:
                with open('possible_coding_insertions.csv', 'a', newline='') as outfile:
                    writer = csv.writer(outfile)
                    writer.writerow([seq, POSITION, gc_content, length, RFP670_seq])




