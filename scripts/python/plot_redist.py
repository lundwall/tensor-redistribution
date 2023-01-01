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
df = pd.concat(df_list)
df = df[df['type'] == 'running'] # remove syncing runs

# Select max across all ranks for each run
df = df.groupby(['id', 'threads', 'chunks', 'mode'])['time'].agg(['max']).reset_index()

# Plot modes for every thread count
for t in range(4, 5):
    plt.figure()
    df_t = df.copy()
    df_t = df_t[df_t['threads'] == t]
    df_t = df_t[df_t['chunks'] == 1]
    fig, axes = plt.subplots()
    sns.violinplot(data = df_t, x = 'mode', y = 'max', ax = axes, inner=None)
    sns.pointplot(data = df_t, x = 'mode', y = 'max', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)
    plt.savefig(f'redist_modes_t{t}.pdf')
    plt.close()

# # Plot threads for manual
plt.figure()
df_m = df.copy()
df_m = df_m[df_m['mode'] == 'manual']
df_m = df_m[df_m['chunks'] == 1]
fig, axes = plt.subplots()
sns.violinplot(data = df_m, x = 'threads', y = 'max', ax = axes, inner=None)
sns.pointplot(data = df_m, x = 'threads', y = 'max', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)
plt.savefig(f'redist_threads_c1.pdf')
plt.close()

# Plot chunks for manual
plt.figure()
df_c = df.copy()
df_c = df_c[df_c['mode'] == 'manual']
df_c = df_c[df_c['threads'] == 4]
fig, axes = plt.subplots()
sns.violinplot(data = df_c, x = 'chunks', y = 'max', ax = axes, inner=None)
sns.pointplot(data = df_c, x = 'chunks', y = 'max', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)
plt.savefig(f'redist_chunks_t4.pdf')
plt.close()