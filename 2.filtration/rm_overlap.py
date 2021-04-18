# Remove all the probes w/ overlaps

# Libraries
import time
import itertools
import numpy as np
import pandas as pd
from datetime import date
from os import listdir
import multiprocessing as mp
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO

def parse_cluster(file_name):
    file = open(file_name, 'r')

    min = None
    index = []
    for l in file:
        if l.startswith('>'):
            fields = l.rstrip().split("|")
            curr_index = int(fields[2])
            index.append(curr_index)
            if not min:
                min = curr_index
            else:
                min = curr_index if curr_index < min else min
        else:
            pass

    min_index = index.index(min)
    index.pop(min_index)
    return index

def parse_dup(file_name):
    file = open(file_name, 'r')

    index = []
    for l in file:
        if l.startswith('>'):
            fields = l.rstrip().split("|")
            curr_index = int(fields[2])
            index.append(curr_index)
        else:
            pass

    return index

if __name__ == '__main__':

    start_time = time.time()

    input_dir = './clust/'
    files = sorted([ f for f in listdir(input_dir) if 'cluster' in f ])
    files_path = map(lambda x: input_dir + x, files)

    table_path = '../../Broad/Data/Probe/pick_all/Marmosets_all_new_probe.xlsx'
    df = pd.read_excel(table_path)

    dup_file = 'dup.fa'
    dup_index = parse_dup(dup_file)

    pool = mp.Pool(processes = (mp.cpu_count()))
    results = pool.map(parse_cluster, files_path)
    pool.close()
    pool.join()

    merged = list(itertools.chain(*results)) + dup_index

    df = df.drop(merged).reset_index(drop = True)

    # Output the header matrix
    excel_name = '../filtered/' + 'Marmosets_filtered_probe.xlsx'

    ## xlsx file
    writer = pd.ExcelWriter(excel_name)

    ## Write DataFrame to a file
    df.to_excel(writer, 'Sheet1', index = False)

    ## Save the result
    writer.save()

    print "--- %s seconds ---" % (time.time() - start_time)
