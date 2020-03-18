#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the translated amino acids.
    """

#Holy shit this was all unneccessary!
#    genetic_code = {'GUC': 'V', 'GUA': 'V', 'GUU': 'V', 'GUG': 'V', 'UAU': 'Y', 'UAC': 'Y', 'UGG': 'W', 
#    'ACG': 'T', 'ACA': 'T', 'ACC': 'T', 'ACU': 'T', 'UCU': 'S', 'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'AGC': 'S', 'AGU': 'S', 
#    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'UUU': 'F', 'UUC': 'F', 'AUG': 'M', 'AAA': 'K', 'AAG': 'K', 
#    'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 
#    'CAU': 'H', 'CAC': 'H', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GGU': 'G', 'GAG': 'E', 'GAA': 'E', 
#    'CAA': 'Q', 'CAG': 'Q', 'UGU': 'C', 'UGC': 'C', 'GAC': 'D', 'GAU': 'D', 'AAC': 'N', 'AAU': 'N',
#    'AGG': 'R', 'AGA': 'R', 'CGG': 'R', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCG': 'A', 'GCA': 'A', 'GCC': 'A', 'GCU': 'A', 
#    'UGA': '*', 'UAG': '*', 'UAA': '*'} 
    
    rna_seq = rna_sequence.upper()
    codon_length = 3
    aa_seq = ""
    stop_codons = ["UGA", "UAG", "UAA"]

    if rna_seq[0:3] not in stop_codons:
        for i in range(0, len(rna_seq), codon_length):
            codon = rna_seq[i:i + 3]
            if codon in stop_codons:
                break
            elif len(codon) != 3:
                break
            else:
                aa_seq += genetic_code[codon]

    return aa_seq


def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    list
        A list of strings; each string is an sequence of amino acids encoded by
        `rna_sequence`.
    """
    seq_list = []
    rna_seq2 = rna_sequence.upper()    
    start_position = 0

    def translate(start_position, rna_seq2, genetic_code):
        codon_length = 3
        aa_seq2 = ""
        stop_codons = ["UGA", "UAG", "UAA"]

        for i in range(start_position, len(rna_seq2), codon_length):
            codon = rna_seq2[i:i + 3]
            if codon in stop_codons:
                break
            elif len(codon) != 3:
                break
            else:
                aa_seq2 += genetic_code[codon]
        print (aa_seq2)

    while start_position < len(rna_seq2):
        start_codon = rna_seq2[start_position:start_position + 3]
        if start_codon == "AUG":
            translated_seq = translate(start_position, rna_seq2, genetic_code)
            seq_list.append(translated_seq)
            start_position += 1
    return seq_list





def get_reverse(sequence):
#    sequence = input()
#    sequence.upper()
    #list(sequence) = 
#    return (sequence[::-1])
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'CGUA'
    >>> get_reverse('ATGC')
    'CGTA'
    """
    #sequence = input()
    seq2 = sequence.upper()
    return(seq2[::-1])

def get_complement(sequence):
    """Get the complement of a `sequence` of nucleotides.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'UACG'
    >>> get_reverse('ATGC')
    'TACG'
    """
    #sequence = input()
    seq3 = sequence.upper()
    trantable=seq3.maketrans("ACTGU", "UGACA")
    return seq3.translate(trantable)

def reverse_and_complement(sequence):
    """Get the reversed and complemented form of a `sequence` of nucleotides.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> reverse_and_complement('AUGC')
    'GCAU'
    >>> reverse_and_complement('ATGC')
    'GCAT'
    """
    seq4a= sequence.upper()
    seq4b = seq4a[::-1]
    
    trantable4=seq4b.maketrans("ACTGU", "UGACA")
    return seq4b.translate(trantable4)

def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """
    pass


if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
