from Bio.pairwise2 import align, format_alignment
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo as matlist

import pathlib
import argparse

def align_seqs(seq1, seq2, seqtype='Protein', seq1name='Sequence_1', seq2name='Sequence_2'):
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # check for errors :)
    if seqtype=='DNA':
        DNA_set = {'A', 'G', 'C', 'T'}
        if not (set(seq1).issubset(DNA_set) and set(seq2).issubset(DNA_set)):
            return 'Your DNA sequences contain unexpected characters. Please ensure they only contain the letters A, G, C, and T.'
    if seqtype=='Protein':
        protein_set = {'I', 'C', 'P', 'Q', 'R', 'E', 'V', 'G', 'F', 'W', 'M', 'H', 'K', 'D', 'S', 'T', 'Y', 'L', 'N', 'A'}
        if not set(seq1).issubset(protein_set) or not set(seq2).issubset(protein_set):
            return 'Your protein sequences contain unexpected characters. Please ensure they only contain single-letter amino acid codes.'
    seq1name = '_'.join(seq1name.split(' '))  # clustal can't have spaces in names
    seq2name = '_'.join(seq2name.split(' '))

    if seqtype == 'Protein':
        matrix = matlist.blosum62
    else:
        matrix = matlist.pam180

    alignment = align.globaldx(seq1, seq2, matrix)[0]
    seqs = {}
    seqs[0] = seq1
    seqs[1] = seq2

    # split into top, symbols, and bottom strings 
    pretty = format_alignment(*alignment, full_sequences=True)
    pretty = pretty.split("\n")

    # find locations of mismatches
    mismatch = [i for i, letter in enumerate(pretty[1]) if letter == '.']

    # find locations of gaps
    gap = [i for i, letter in enumerate(pretty[1]) if letter == ' ']

    # replace spaces with the no blank space whatever it's called char 
    # so that spaces show in the html rendering 
    symbols_list = list(pretty[1])
    for i in range(len(symbols_list)):
        if symbols_list[i] == ' ':
            symbols_list[i] = '&nbsp;'

    # write a string in the clustalw file format (https://meme-suite.org/meme/doc/clustalw-format.html)
    clustal = 'CLUSTAL W\n\n'

    # calculate spacing between names and start of actual sequences' alignment info 
    max_name_len = max(len(seq1name), len(seq2name))
    spacers = [max_name_len - len(seq1name) + 6, 
               max_name_len - len(seq2name) + 6, 
               max_name_len + 6]
    
    i = 1
    while 60*i <= len(pretty[1]):
        clustal += seq1name + ' ' * spacers[0] + pretty[0][60*(i-1):60*i] + f' {60*i}' + '\n' + \
                   seq2name + ' ' * spacers[1] + pretty[2][60*(i-1):60*i] + f' {60*i}' + '\n' + \
                   ' ' * spacers[2] + pretty[1].replace('|', '*')[60*(i-1):60*i] + '\n\n'
        i += 1
    
    # the last bit
    clustal += seq1name + ' ' * spacers[0] + pretty[0][60*(i-1):] + f' {len(pretty[1])}' + '\n' + \
                   seq2name + ' ' * spacers[1] + pretty[2][60*(i-1):] + f' {len(pretty[1])}' + '\n' + \
                   ' ' * spacers[2] + pretty[1].replace('|', '*')[60*(i-1):] + '\n\n'

    # determine indices of bases/AAs, excluding gaps!!
    top_index, bottom_index = ['&nbsp;'] * len(symbols_list), ['&nbsp;'] * len(symbols_list)

    counter10 = 0
    for i, char in enumerate(list(pretty[0])):
        if char != '-':
            counter10 += 1
            if counter10 % 10 == 0:
                top_index[i] = str(counter10)

    
    counter10 = 0
    for i, char in enumerate(list(pretty[2])):
        if char != '-':
            counter10 += 1
            if counter10 % 10 == 0:
                bottom_index[i] = str(counter10)


    for i in range(len(top_index)-1, 0, -1):
        char = top_index[i]
        if char != '&nbsp;':
            num_digits = len(char)   # remove spaces after double/triple/ digit numbers for formatting
            for _ in range(num_digits - 1):
                try: # account for edge case where there's nothing after the number
                    del top_index[i+1]  
                except:
                    pass

    for i in range(len(bottom_index)-1, 0, -1):
        char = bottom_index[i]
        if char != '&nbsp;':
            num_digits = len(char)   # remove spaces after double/triple/ digit numbers for formatting
            for _ in range(num_digits - 1):
                try: # account for edge case where there's nothing after the number
                    del bottom_index[i+1]  
                except:
                    pass

    output = {}
    output['alignment-top'] = [seq1name] + ['&nbsp;'] * spacers[0] + list(pretty[0])
    output['alignment-symbols'] = ['&nbsp;'] * spacers[2] + symbols_list
    output['alignment-bottom'] = [seq2name] + ['&nbsp;'] * spacers[1] + list(pretty[2])
    output['top-index'] = ''.join(['&nbsp;'] * spacers[2] + top_index)
    output['bottom-index'] = ''.join(['&nbsp;'] * spacers[2] + bottom_index)
    output['mismatch'] = [i + spacers[2] for i in mismatch]  # seems like these are more conservative mismatches?
    output['gap'] = [i + spacers[2] for i in gap]   # gaps and less conservative mismatches?
    output['clustal'] = clustal
    output['spacer'] = spacers[2]  # hella goofy ik but necessary for my silly lil code

    return output

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("seqfile1")
    parser.add_argument("seqfile2")

    args = parser.parse_args()

    align(args.seqfile1, args.seqfile2)