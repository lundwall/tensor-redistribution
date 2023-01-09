import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import sys
import os

folder = sys.argv[1]
size_to_df = {}
for filename in os.listdir(os.getcwd() + '/' + folder):
    df = pd.read_csv(folder + '/' + filename, delim_whitespace=True, comment='#')
    if 'chunk' not in df:
        append_column = [1 for x in range(len(df.index))]
        df = df.assign(chunk=append_column)
    name_elements = filename.split('.')[1].split('_')
    size_elements = name_elements[-5::]
    try:
        size_tuple = (int(size_elements[0]), int(size_elements[1]), int(size_elements[2]), int(size_elements[3]), int(size_elements[4]))
        if size_tuple in size_to_df:
            size_to_df[size_tuple].append(df)
        else:
            size_to_df[size_tuple] = [df]
    except ValueError:
        print("some size is not integer")

total_df = []
for key, value in size_to_df.items():
    size_to_df[key] = pd.concat(value)
    size_string = str(key[0]) + "*" + str(key[1]) + "*" + str(key[2]) + "*" + str(key[3]) + "*" + str(key[4])
    append_column = [size_string for x in range(len(size_to_df[key].index))]
    size_to_df[key] = size_to_df[key].assign(size_string=append_column)
    total_size = key[0] * key[1] * key[2] * key[3] * key[4] * 4 / 1000
    append_column = [total_size for x in range(len(size_to_df[key].index))]
    size_to_df[key] = size_to_df[key].assign(total_size=append_column)
    total_df.append(size_to_df[key])

df_all = pd.concat(total_df)
df_all = df_all[df_all['type'] == 'running'] # remove syncing runs

# manual packing
df_t1_c1 = df_all[(df_all['chunk'] == 1) & (df_all['threads'] == 1)]
df_t1_c1 = df_t1_c1[(df_t1_c1['mode'] == 'datatype') | (df_t1_c1['mode'] == 'manual')]
df_t1_c1 = df_t1_c1.groupby(['id', 'mode', 'size_string', 'total_size'])['time'].agg(['max']).reset_index()
df_t1_c1 = df_t1_c1.sort_values(['total_size'])

plt.figure()
fig, axes = plt.subplots()
sns.violinplot(data = df_t1_c1, x = 'size_string', y = 'max', hue = 'mode', ax = axes)
# sns.pointplot(data = df_t1_c1, x = 'mode', y = 'max', order=sorted(df['mode'].unique()), ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)

plt.savefig(f'manual_vs_datatype.pdf')
plt.close()

# manual packing with thread acceleration
df_c1_lsize = df_all[((df_all['size_string'] == '24*24*24*24*24') & (df_all['mode'] == 'manual') & (df_all['chunk'] == 1)) | ((df_all['size_string'] == '24*24*24*24*24') & (df_all['mode'] == 'datatype') & (df_all['threads'] == 1))]
df_c1_lsize = df_c1_lsize.groupby(['id', 'mode', 'threads'])['time'].agg(['max']).reset_index()
plt.figure()
fig, axes = plt.subplots()
sns.violinplot(data = df_c1_lsize, x = 'threads', y = 'max', order=sorted(df_c1_lsize['threads'].unique()), hue = 'mode',ax = axes)
# sns.pointplot(data = df_c1_lsize, x = 'chunk', y = 'max', order=sorted(df_c1_lsize['threads'].unique()), hue = 'mode', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)

plt.savefig(f'manual_with_threads_vs_datatype.pdf')
plt.close()

# manual packing with chunking
df_t1_lsize = df_all[(df_all['size_string'] == '24*24*24*24*24') & (df_all['threads'] == 1)]
df_t1_lsize = df_t1_lsize.groupby(['id', 'mode', 'chunk'])['time'].agg(['max']).reset_index()
df_t1_lsize = df_t1_lsize[(df_t1_lsize['mode'] == 'datatype') | (df_t1_lsize['mode'] == 'manual')]
plt.figure()
fig, axes = plt.subplots()
sns.violinplot(data = df_t1_lsize, x = 'chunk', y = 'max', hue = 'mode', ax = axes)
# sns.pointplot(data = df_t1_lsize, x = 'chunk', y = 'max', order=sorted(df_t1_lsize['chunk'].unique()), hue = 'mode', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)

plt.savefig(f'manual_packing_with_chunking_vs_datatype.pdf')
plt.close()

# MPI one sided
df_c1_lsize = df_all[(df_all['size_string'] == '24*24*24*24*24') & (df_all['chunk'] == 1) & (df_all['threads'] == 1)]
df_c1_lsize = df_c1_lsize.groupby(['id', 'mode'])['time'].agg(['max']).reset_index()
plt.figure()
fig, axes = plt.subplots()
sns.violinplot(data = df_c1_lsize, x = 'mode', y = 'max', ax = axes)
plt.savefig(f'manual_packing_with_one_sided_vs_datatype.pdf')
plt.close()