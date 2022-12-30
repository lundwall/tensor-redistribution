import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import sys
import os

df_list = []
columns = ['mode', 'rank', 'time', 'type', 'threads']
folder = sys.argv[1]
for filename in os.listdir(os.getcwd() + '/' + folder):
    df = pd.read_csv(folder + '/' + filename, delim_whitespace=True, comment='#', usecols=columns)
    df = df[df['type'] == 'running']
    df_list.append(df)
df = pd.concat(df_list)

df_manual = df[df['mode'] == 'manual']
df_datatype = df[df['mode'] == 'datatype']

max_rank_manual = df_manual.groupby(['rank'])['time'].median().idxmax()
max_rank_datatype = df_datatype.groupby(['rank'])['time'].median().idxmax()

df_manual = df_manual[df_manual['rank'] == max_rank_manual]
df_datatype = df_datatype[df_datatype['rank'] == max_rank_manual]

df = pd.concat([df_manual, df_datatype])

fig, axes = plt.subplots()
sns.violinplot(data = df, x = 'mode', y = 'time', ax = axes, inner=None)
sns.pointplot(data = df, x = 'mode', y = 'time', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)

plt.savefig('density_t4_direct_redist.pdf')