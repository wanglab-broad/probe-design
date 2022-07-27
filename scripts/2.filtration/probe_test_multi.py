#!/usr/bin/python

# Functions
import time
import melting
import numpy as np
import pandas as pd
from datetime import date
from os import listdir
import multiprocessing as mp
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO

def ReverseComplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def split_primer(seq, offset=4):
    """ Splits primer :seq: into 2 of equal Tm. Search for +/- :offset: around midpoint. """
    N = len(seq)
    midpoint = int(round(float(N)/2))
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
        if len(front_seq) >= 19 and len(front_seq) <= 25 and len(back_seq) >= 19 and len(back_seq) <= 25:
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

    return f,b, temp_pairs[min_idx], temp_diffs[min_idx], spacer#, len(f), len(b)

def parse_picky_new(x):
    fname, source_index = x
    columns = ['Gene', 'Id', 'Start', 'End', 'TargetSequence', 'Tm', 'ReverseComplement', 'Primer', 'Primer Tm', 'Spacer', 'Padlock', 'Padlock Tm', 'DeltaTm']
    df = pd.DataFrame(columns=columns)

    f = open(fname,'r')
    k = 0 # for each record, identify line
    gene_seq = ""
    curr_df = pd.DataFrame(columns=columns)
    curr_line = 0 # current line number of the output dataframe
    has_shared = False

    # Iterate through lines in pick file
    for l in f:
        fields = l.rstrip().split("\t") # remove trailing space and split by tab
        # print fields
        k += 1

        # Check content in each line
        if fields == ['']: # empty line
            if not has_shared:
                df = df.append(curr_df, ignore_index = True)

            k = 0
            curr_df = pd.DataFrame(columns=columns)
            curr_line = 0
            has_shared = False

        elif k == 1: # ['ACTGGGATGTTCGGAGCATTCAACGCTGGTTTCGACAAAGAC', '42', '1', '1', '83.51', '52.87']
            gene_seq = fields[0]
            primer,padlock,Tm_pairs,deltaT,spacer = split_primer(ReverseComplement(gene_seq))
        elif k > 1 and fields[0] == "U": # ['U', '83.51', '0', '41', '319', '360', '16086.1|Cers6']
            Tm = fields[1]
            start = int(fields[4])
            stop = int(fields[5])
            (seq_id, name) = fields[6].split("|")
            curr_row = [name, seq_id, start, stop, gene_seq, Tm, ReverseComplement(gene_seq), primer, Tm_pairs[0], spacer, padlock, Tm_pairs[1], deltaT]
            #print curr_row
            curr_df.loc[curr_line] = curr_row
            curr_line += 1
        elif k > 1 and fields[0] == "S": # ['S', '84.74', '0', '39', '240', '279', 'CCDS10674.1|SULT1A3']
            has_shared = True
        elif fields[0] == "<" or fields[0] == ">":
            k = -1

    df = df.sort_values(['Gene', 'DeltaTm']).reset_index(drop = True) ## Sort
    df['source_index'] = [source_index] * len(df.index)
    f.close()

    return df

def make_seq_record(seq, name):
    seq_record = SeqRecord(Seq(seq), id = name, description = '', letter_annotations = {"phred_quality":[40] * len(seq)})
    return seq_record

def output_log(log_file, total_probe, dup_record, pre_sns, sns_percent, final_probe):
    print >> log_file, """\n
Total number of probes from PICKY files: %d
Number of duplicated records removed: %d
Before removing sns, there are %d probes.
%d%% probes w/ sns have been removed.
Total number of remainning probes after pre-processing : %d \n""" % (total_probe, dup_record, pre_sns, sns_percent, final_probe)



if __name__ == '__main__':
    # Iterate through species
    for s in ['Marmosets']:
        species = s
        start_time = time.time()
        print species
        # Parameter
        # species = 'H_sapiens'
        # species = 'M_musculus'
        # species = 'Marmosets'

        # Load Picky2.0 output
        input_path = 'Broad/Data/Picky/' + species + '/new'
        log_name = 'res/' + species + '_log.txt'
        log_file = open(log_name, 'w')

        files = sorted([ f for f in listdir(input_path) if '.picky' in f])
        # print files

        # Construct the probe picking order:
        # #1 ccds #1 rna
        # #2 ccds #2 rna
        # #3 ccds #3 rna
        # will choose the probe from ccds/more strict set
        source_index = 1
        p = []

        # Marmosets dataset doesn't have ccds
        if species == 'Marmosets':
            order = len(files)
            ccds = False
        else:
            order = len(files)/2
            ccds = True

        for i in range(order):
            curr_ccds = input_path + '/' + files[i]
            p.append((curr_ccds, source_index))
            source_index += 1

            if ccds:
                curr_rna = input_path + '/' + files[i + 3]
                p.append((curr_rna, source_index))
                source_index += 1

        # Construct multi-processing pool
        print '== File == Index =='
        for process in p:
            print process[0] + '  ' + str(process[1])
            print >> log_file, process[0] + '  ' + str(process[1])

        pool = mp.Pool(processes = (mp.cpu_count()))
        results = pool.map(parse_picky_new, p)
        pool.close()
        pool.join()

        df = pd.concat(results)

        # print results_df.shape
        df = df.sort_values(['Gene', 'source_index'])

        total_probe = len(df.index)

        colnames = df.columns.values.tolist()
        #print colnames
        df_dp = df.drop_duplicates(colnames[:len(colnames) - 1])
        df_dp = df_dp.drop_duplicates(['Gene', 'TargetSequence']) # remove identical probes within gene group
        df_dp = df_dp.drop_duplicates(['TargetSequence'], keep = False) # remove duplicated probes

        dup_record = len(df.index) - len(df_dp.index)

        old_nrow = len(df_dp.index)
        sns_drop = ['AAAAA', 'TTTTT', 'GGGGG', 'CCCCC']
        df_dp = df_dp[~df_dp['ReverseComplement'].str.contains('|'.join(sns_drop))]
        af_nrow = len(df_dp.index)
        af_percent = round((old_nrow - af_nrow) / float(old_nrow), 2) * 100
        print 'There are %d%% probes filtered out because of continuous single nucleotide sequence!' % af_percent
        # df_dp = df_dp.groupby('Gene').head(3).reset_index(drop = True)

        df_dp = df_dp.reset_index(drop = True)
        df_dp['temp_index'] = df_dp.index.astype(str)
        df_dp['fasta_id'] = df_dp[['Gene', 'Id', 'temp_index']].apply(lambda x: '|'.join(x), axis=1)
        df_dp = df_dp.drop('temp_index', 1)

        final_probe = len(df_dp.index)
        # print df_dp
        # print df_dp[['ReverseComplement', 'Primer', 'Spacer', 'Padlock']]

        seq_records = map(make_seq_record, df_dp.TargetSequence, df_dp.fasta_id)

        # for index, row in df_dp.iterrows():
        #     print row['Id'], row['Start'], row['End']

        # Output FASTA
        fasta_name = 'res/' + species + '_all_new_probe.fa'
        out_file = open(fasta_name, 'w')
        Bio.SeqIO.write(seq_records, out_file, 'fasta')

        # Output FASTQ w/ fake quality score
        fastq_name = 'res/' + species + '_all_new_probe.fq'
        out_file = open(fastq_name, 'w')
        Bio.SeqIO.write(seq_records, out_file, 'fastq')

        # Output the header matrix
        excel_name = 'res/' + species + '_all_new_probe.xlsx'

        ## xlsx file
        writer = pd.ExcelWriter(excel_name)

        ## Write DataFrame to a file
        df_dp.to_excel(writer, 'Sheet1', index = False)

        ## Save the result
        writer.save()

        # Output log
        output_log(log_file, total_probe, dup_record, old_nrow, af_percent, final_probe)

        print "--- %s seconds ---" % (time.time() - start_time)
        print >>log_file, "--- %s seconds ---" % (time.time() - start_time)
