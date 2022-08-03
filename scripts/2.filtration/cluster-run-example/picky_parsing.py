#!/usr/bin/env python
# coding: utf-8

# cml input
import sys, os

# get input files
input_path = sys.argv[1]
files = sorted([ f for f in os.listdir(input_path) if '.picky' in f])
print(files)

# import probe design module
from probe import *
import timeit

# start the timer
start = timeit.default_timer()

# setup multiple proessing for parsing the files
p = [] # (input_file_path, source_index)
for i, current_file in enumerate(files):
    current_tuple = (os.path.join(input_path, current_file), i+1)
    p.append(current_tuple)

pool = mp.Pool(processes = (mp.cpu_count()))
results = pool.starmap(parse_picky_new, p)
pool.close()
pool.join()

df_combined = pd.concat(results)

# set output path
output_path = os.path.join(input_path, 'output')
if not os.path.exists(output_path):
    os.mkdir(output_path)

# save results
df_combined.to_csv(os.path.join(output_path, 'parsed_results.csv'))

# save log (csv)
stop = timeit.default_timer()
run_time = round((stop - start) / 60, 2)

with open(os.path.join(output_path, 'log.txt'), 'w') as log_file:
    print(f"Number of records: {df_combined.shape[0]}", file=log_file)
    print(f"Run time: {run_time} min", file=log_file)
