import pandas as pd
from os import listdir

input_dir = '../../Broad/Data/Probe/filtered/'
files = sorted([ f for f in listdir(input_dir) if 'filtered' in f ])
files_path = map(lambda x: input_dir + x, files)

for p in range(len(files_path)):
    df = pd.read_excel(files_path[p])
    df = df.groupby('Gene').head(3).reset_index(drop = True)

    curr_file = files[p].split('.')
    excel_name = input_dir + curr_file[0] + '_head3.xlsx'
    print excel_name
    writer = pd.ExcelWriter(excel_name)
    df.to_excel(writer, 'Sheet1', index = False)
    writer.save()
