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
df = df.groupby(['id', 'threads', 'chunks', 'size', 'mode'])['time'].agg(['max']).reset_index()

# Plot modes for t1, c1 for each size
for s in sorted(df['size'].unique()):
    plt.figure()
    df_t = df.copy()
    df_t = df_t[df_t['threads'] == 1]
    df_t = df_t[df_t['chunks'] == 1]
    df_t = df_t[df_t['size'] == s]
    fig, axes = plt.subplots()
    sns.violinplot(data = df_t, x = 'mode', y = 'max', ax = axes, inner=None)
    sns.pointplot(data = df_t, x = 'mode', y = 'max', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)
    plt.savefig(f'redist_modes_t1_c1_s{int(s / 1e6)}MB.pdf')
    plt.close()

# Plot threads for manual for every chunk
for c in sorted(df['chunks'].unique()):
    plt.figure()
    df_m = df.copy()
    df_m = df_m[df_m['mode'] == 'manual']
    df_m = df_m[df_m['chunks'] == c]
    df_m = df_m[df_m['size'] == 234375000]
    fig, axes = plt.subplots()
    sns.violinplot(data = df_m, x = 'threads', y = 'max', ax = axes, inner=None)
    sns.pointplot(data = df_m, x = 'threads', y = 'max', ax = axes, estimator=np.median, color='black', join=True, scale=0.5, errwidth=1, capsize=0.01)
    plt.savefig(f'redist_threads_c{c}.pdf')
    plt.close()

# Plot chunks for manual for every thread count
for t in sorted(df['threads'].unique()):
    plt.figure()
    df_c = df.copy()
    df_c = df_c[df_c['mode'] == 'manual']
    df_c = df_c[df_c['threads'] == t]
    df_c = df_c[df_c['size'] == 234375000]
    fig, axes = plt.subplots()
    # sns.violinplot(data = df_c, x = 'chunks', y = 'max', ax = axes, inner=None)
    sns.pointplot(data = df_c, x = 'chunks', y = 'max', ax = axes, estimator=np.median, color='black', join=True, scale=0.5, errwidth=1, capsize=0.01)
    plt.savefig(f'redist_chunks_t{t}.pdf')
    plt.close()

# Chunks/Threads heatmap for each size
for s in sorted(df['size'].unique()):
    plt.figure()
    df_h = df.copy()
    df_h = df_h[df_h['mode'] == 'manual']
    df_h = df_h[df_h['size'] == s]
    # Flatten across ids by taking median
    df_h = df_h.groupby(['threads', 'chunks'])['max'].agg(['median']).reset_index()
    df_h = df_h.pivot(index='chunks', columns='threads', values='median')
    fig, axes = plt.subplots()
    sns.heatmap(df_h)
    plt.savefig(f'redist_heatmap_s{int(s / 1e6)}MB.pdf')
    plt.close()