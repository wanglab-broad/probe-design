# The python codebase for the probe design in Wang Lab
# Jiahao Huang - 2022

# libraries
import os
import time
import numpy as np
import pandas as pd
import melting
import multiprocessing as mp
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO


# Get first n records from a dataframe based on specific field
def get_first_n_records(df: pd.DataFrame, field: str, n: int) -> pd.DataFrame:
    """

    Parameters
    ----------
    df: input pandas dataframe
    field: the column used to group the dataframe
    n: number of records needed

    Returns
    -------
    df: output pandas dataframe

    """
    return df.groupby(field).head(n).reset_index(drop=True)


# Get the reverse complement of an input sequence
def reverse_complement(seq):
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return "".join([seq_dict[base] for base in reversed(seq)])


# Split the sequence into two pieces
def split_primer(seq, offset=4):
    """ Splits primer :seq: into 2 of equal Tm. Search for +/- :offset: around midpoint. """
    n = len(seq)
    midpoint = int(round(float(n)/2))
    front_seqs = []
    back_seqs = []
    temp_diffs = []
    temp_pairs = []
    # spacers = []
    max_offset_plus = offset
    max_offset_minus = offset
    offsets = range(-offset, offset)
    for i in offsets:
        front_seq = seq[:(midpoint+i)]
        back_seq = seq[(midpoint+i):]
        if len(seq) == 40 or len(seq) == 41:
            # spacers.append(front_seq[len(front_seq)-1])
            front_seq = front_seq[:(len(front_seq)-1)]
        else:
            # spacers.append(front_seq[len(front_seq)-1] + back_seq[0])
            front_seq = front_seq[:(len(front_seq)-1)]
            back_seq = back_seq[1:]
        front_tm = melting.temp(front_seq, Na_c=300, Mg_c=0)
        back_tm = melting.temp(back_seq, Na_c=300, Mg_c=0)
        if 19 <= len(front_seq) <= 25 and 19 <= len(back_seq) <= 25:
            front_seqs.append(front_seq)
            back_seqs.append(back_seq)
            temp_diffs.append(abs(front_tm-back_tm))
            temp_pairs.append((front_tm, back_tm))
    temp_diffs = np.array(temp_diffs)
    min_idx = np.where(temp_diffs == np.min(temp_diffs))[0][0]
    f = front_seqs[min_idx]
    b = back_seqs[min_idx]

    sp_loc = seq.find(b)
    if len(seq) == 40 or len(seq) == 41:
        spacer = seq[sp_loc-1:sp_loc]
    else:
        spacer = seq[sp_loc-2:sp_loc]

    return f, b, temp_pairs[min_idx], temp_diffs[min_idx], spacer #, len(f), len(b)


# Parse single picky output (.picky) file
def parse_picky_new(fname, source_index):
    """
    Parameters
    ----------
    fname: input path of the picky file
    source_index: an integer assigned to the current picky file

    Returns
    -------
    df: output pandas dataframe with parsed picky file

    """
    # Set the columns for the dataframe output
    columns = ['Annotation',
               'Start', 'End', 'TargetSequence',
               'Tm', 'ReverseComplement',
               'Primer', 'Primer Tm', 'Spacer',
               'Padlock', 'Padlock Tm', 'DeltaTm']

    df = pd.DataFrame(columns=columns)

    # Read the file
    f = open(fname, 'r')
    k = 0  # for each record, identify line
    gene_seq = ""
    current_df = pd.DataFrame(columns=columns)
    current_line = 0  # current line number of the output dataframe
    has_shared = False

    # Iterate through lines in pick file
    for l in f:
        fields = l.rstrip().split("\t") # remove trailing space and split by tab
        # print fields
        k += 1

        # Check content in each line
        if fields == ['']:  # empty line
            if not has_shared:
                df = pd.concat([df, current_df], ignore_index=True)
                # df = df.append(current_df, ignore_index = True)

            # Reset parameters
            k = 0
            current_df = pd.DataFrame(columns=columns)
            current_line = 0
            has_shared = False

        elif k == 1:  # ['ACTGGGATGTTCGGAGCATTCAACGCTGGTTTCGACAAAGAC', '42', '1', '1', '83.51', '52.87']
            gene_seq = fields[0]
            primer, padlock, Tm_pairs, deltaT, spacer = split_primer(reverse_complement(gene_seq))

        elif k > 1 and fields[0] == "U":  # ['U', '83.51', '0', '41', '319', '360', '16086.1|Cers6']
            Tm = fields[1]
            start = int(fields[4])
            stop = int(fields[5])
            annot = fields[6]
            current_row = [annot,
                           start, stop, gene_seq,
                           Tm, reverse_complement(gene_seq),
                           primer, Tm_pairs[0], spacer,
                           padlock, Tm_pairs[1], deltaT]
            current_df.loc[current_line] = current_row
            current_line += 1

        elif k > 1 and fields[0] == "S":  # ['S', '84.74', '0', '39', '240', '279', 'CCDS10674.1|SULT1A3']
            has_shared = True

        elif fields[0] == "<" or fields[0] == ">":
            k = -1

    df = df.sort_values(['Annotation', 'DeltaTm']).reset_index(drop=True) ## Sort
    df['source_index'] = [source_index] * len(df.index)
    f.close()

    return df


# Generate sequence record
def make_seq_record(seq, name):
    seq_record = SeqRecord(Seq(seq), id=name, description='', letter_annotations={"phred_quality": [40] * len(seq)})
    return seq_record
