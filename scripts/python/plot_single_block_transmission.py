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
        size_string = size_elements[0] + "_" + size_elements[1] + "_" + size_elements[2] + "_" + size_elements[3] + "_" + size_elements[4]
        if size_string in size_to_df:
            size_to_df[size_string].append(df)
        else:
            size_to_df[size_string] = [df]
    except ValueError:
        print("some size is not integer")

total_df = []
for key, value in size_to_df.items():
    size_to_df[key] = pd.concat(value)
    append_column = [key for x in range(len(size_to_df[key].index))]
    size_to_df[key] = size_to_df[key].assign(size=append_column)
    total_df.append(size_to_df[key])

df_all = pd.concat(total_df)
df_all = df_all[df_all['type'] == 'running'] # remove syncing runs

df_t1_c1 = df_all[((df_all['chunk'] == 1) & (df_all['threads'] == 1)) | ((df_all['chunk'].isnull()) & (df_all['threads'] == 1))]
df_t1_c1 = df_t1_c1.groupby(['id', 'mode', 'size'])['time'].agg(['max']).reset_index()

# For size, chunk = 1, thread = 1
for key, value in size_to_df.items():
    df = df_t1_c1[df_t1_c1['size'] == key]
    plt.figure()
    fig, axes = plt.subplots()
    sns.violinplot(data = df, x = 'mode', y = 'max', order=sorted(df['mode'].unique()), ax = axes)
    # sns.pointplot(data = df, x = 'mode', y = 'max', order=sorted(df['mode'].unique()), ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)

    plt.savefig(f'size_{key}.pdf')
    plt.close()

# Largest size, plot chunk, thread = 1
df_t1_lsize = df_all[(df_all['size'] == '24_24_24_24_24') & (df_all['threads'] == 1)]
df_t1_lsize = df_t1_lsize.groupby(['id', 'mode', 'chunk'])['time'].agg(['max']).reset_index()
plt.figure()
fig, axes = plt.subplots()
sns.violinplot(data = df_t1_lsize, x = 'chunk', y = 'max', order=sorted(df_t1_lsize['chunk'].unique()), hue = 'mode',ax = axes)
# sns.pointplot(data = df_t1_lsize, x = 'chunk', y = 'max', order=sorted(df_t1_lsize['chunk'].unique()), hue = 'mode', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)

plt.savefig(f'lsize_chunk_24_24_24_24_24.pdf')
plt.close()

# Largest size2, plot chunk, thread = 1
df_t1_lsize = df_all[(df_all['size'] == '12_12_12_12_96') & (df_all['threads'] == 1)]
df_t1_lsize = df_t1_lsize.groupby(['id', 'mode', 'chunk'])['time'].agg(['max']).reset_index()
plt.figure()
fig, axes = plt.subplots()
sns.violinplot(data = df_t1_lsize, x = 'chunk', y = 'max', order=sorted(df_t1_lsize['chunk'].unique()), hue = 'mode',ax = axes)
# sns.pointplot(data = df_t1_lsize, x = 'chunk', y = 'max', order=sorted(df_t1_lsize['chunk'].unique()), hue = 'mode', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)

plt.savefig(f'lsize_chunk_12_12_12_12_96.pdf')
plt.close()

# Largest size, plot thread, chunk = 1
df_c1_lsize = df_all[(df_all['size'] == '24_24_24_24_24') & (df_all['chunk'] == 1)]
df_c1_lsize = df_c1_lsize.groupby(['id', 'mode', 'threads'])['time'].agg(['max']).reset_index()
plt.figure()
fig, axes = plt.subplots()
sns.violinplot(data = df_c1_lsize, x = 'threads', y = 'max', order=sorted(df_c1_lsize['threads'].unique()), hue = 'mode',ax = axes)
# sns.pointplot(data = df_c1_lsize, x = 'chunk', y = 'max', order=sorted(df_c1_lsize['threads'].unique()), hue = 'mode', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)

plt.savefig(f'lsize_threads_24_24_24_24_24.pdf')
plt.close()

# Largest size2, plot chunk, chunk = 1
df_c1_lsize = df_all[(df_all['size'] == '12_12_12_12_96') & (df_all['chunk'] == 1)]
df_c1_lsize = df_c1_lsize.groupby(['id', 'mode', 'threads'])['time'].agg(['max']).reset_index()
plt.figure()
fig, axes = plt.subplots()
sns.violinplot(data = df_c1_lsize, x = 'threads', y = 'max', order=sorted(df_c1_lsize['threads'].unique()), hue = 'mode',ax = axes)
# sns.pointplot(data = df_c1_lsize, x = 'chunk', y = 'max', order=sorted(df_c1_lsize['threads'].unique()), hue = 'mode', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)

plt.savefig(f'lsize_threads_12_12_12_12_96.pdf')
plt.close()