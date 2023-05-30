#!/usr/bin/env python3
# Import pysam module
import pysam

# Open the SAM file
samfile = pysam.AlignmentFile("-", "rb")

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return ''.join([complement[base] for base in dna[::-1]])

# Loop over each read in the file
for read in samfile:
    #df = pd.DataFrame()
    # Get the ML tag value
    try:
        ml_tag = read.get_tag("ML")
    except:
        continue

    # If it has, print the ML tag value and the read name
    if ml_tag:
        # Get the read sequence
        seq = read.query_sequence
    
        # If the read is mapped to the reverse strand, reverse complement the sequence
        if read.is_reverse:
            seq = reverse_complement(seq)
    
        # Find all the positions of CG motifs in the sequence
        cg_pos = []
        for i in range(len(seq) - 1):
            if seq[i:i+2] == "CG":
                cg_pos.append(i)
        # Print the ML tag value and the positions of CG motifs
        read.query_name
        LIST=zip(cg_pos, ml_tag)
        for i,j in LIST:
            print(f"{read.query_name}\t{i}\t{j}")
