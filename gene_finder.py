# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Sampei Omichi

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # Complete
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == "G":
        return 'C'

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # Complete
    reverse_complement = ''
    for i in dna:
        reverse_complement = get_complement(i) + reverse_complement
    return reverse_complement

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # Complete
    stop_codon = ['TAG', 'TAA', 'TGA']
    dna = dna[3:]
    for i in range(len(dna)):
        if dna[i:i+3] in stop_codon:
            return 'ATG' + dna[0:i]

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # Complete
    stop_codon_list = ['TAG', 'TAA', 'TGA']
    ORF = ''
    all_ORFs_oneframe = []
    start_codon = None
    stop_codon = None
    for i in range(len(dna)):
        if i % 3 == 0 and dna[i:i+3] == 'ATG':
            if start_codon != None and stop_codon == None:
                pass
            else:
                start_codon = i
        if i % 3 == 0 and dna[i:i+3] in stop_codon_list:
            stop_codon = i
        if i % 3 == 0 and dna[i:i+3] and start_codon != None and stop_codon != None:
            ORF = dna[start_codon: stop_codon]
        if ORF != '':
            all_ORFs_oneframe.append(ORF)
            start_codon = None
            stop_codon = None
            ORF = ''
        elif i == len(dna)-1 and start_codon != None:
            ORF = dna[start_codon:i+1]
            all_ORFs_oneframe.append(ORF)
    return all_ORFs_oneframe

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # Complete
    all_ORFS = []
    frame_1 = find_all_ORFs_oneframe(dna)
    frame_2 = find_all_ORFs_oneframe(dna[1:])
    frame_3 = find_all_ORFs_oneframe(dna[2:])
    ORFS =  frame_1 + frame_2 + frame_3
    for i in ORFS:
        if i not in all_ORFS and i[0:3] == 'ATG':
            all_ORFS.append(i)
    return all_ORFS

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # Complete
    front = dna
    reverse = get_reverse_complement(dna)
    ORFs_both_strands = find_all_ORFs(front) + find_all_ORFs(reverse)
    return ORFs_both_strands

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # Complete
    ORFs_both_strands = find_all_ORFs_both_strands(dna)
    if len(ORFs_both_strands) > 0:
        return max(ORFs_both_strands, key=len)
    else:
        pass

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # Complete
    longest_ORFs = []
    random_dna = ''
    for i in range(num_trials):
        random_dna = shuffle_string(dna)
        ORF = longest_ORF(random_dna)
        if ORF != None:
            longest_ORFs.append(ORF)
    if len(longest_ORFs) > 0:
        return max(longest_ORFs, key=len)

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # Complete
    amino_acids = ''
    while dna:
        triple_codons = dna[:3]
        if len(triple_codons) == 3:
            amino_acids = amino_acids + aa_table[triple_codons]
        dna = dna[3:]
    return amino_acids

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # Complete
    threshold = longest_ORF_noncoding(dna, 1500)
    all_ORFs = find_all_ORFs_both_strands(dna)
    amino_acids_list = []
    for i in all_ORFs:
        if len(i) < len(threshold):
            del all_ORFs[all_ORFs.index(i)]
    for i in all_ORFs:
        amino_acids_list.append(coding_strand_to_AA(i))
    return amino_acids_list

dna = load_seq("./data/X73525.fa")
print(gene_finder(dna))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
