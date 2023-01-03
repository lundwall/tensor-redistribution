import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import sys
import os

df_list = []
folder = sys.argv[1]
for filename in os.listdir(os.getcwd() + '/' + folder):
    df = pd.read_csv(folder + '/' + filename, delim_whitespace=True, comment='#')
    df_list.append(df)
df_saved = pd.concat(df_list)
df_saved = df_saved[df_saved['type'] == 'running'] # remove syncing runs

# Select max across all ranks for each run
df = df_saved.groupby(['id', 'threads', 'mode'])['time'].agg(['max']).reset_index()

for t in sorted(df['threads'].unique()):
    plt.figure()
    df = df_saved.copy()
    print(df)
    df = df[df['threads'] == t] # select threads from loop
    # df = df[df['mode'] != 'template'] # template and put_datatype are too big to be nice, so they're excluded
    # df = df[df['mode'] != 'put_datatype']

    fig, axes = plt.subplots()
    sns.violinplot(data = df, x = 'mode', y = 'time', order=sorted(df['mode'].unique()), ax = axes, inner=None)
    sns.pointplot(data = df, x = 'mode', y = 'time', order=sorted(df['mode'].unique()), ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)

    plt.savefig(f'single_t{t}.pdf')
    plt.close()