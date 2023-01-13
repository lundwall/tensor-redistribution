import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib import ticker as mticker
import numpy as np
import sys
import os

file_suffix = '_final'

df_list = []
folder = sys.argv[1]
for filename in os.listdir(os.getcwd() + '/' + folder):
    df = pd.read_csv(folder + '/' + filename, delim_whitespace=True, comment='#')
    df_list.append(df)
df = pd.concat(df_list)
df = df[df['type'] == 'running'] # remove syncing runs

# Select max across all ranks for each run
df = df.groupby(['id', 'threads', 'chunks', 'size', 'mode'])['time'].agg(['max']).reset_index()

print(df.groupby(['threads', 'chunks', 'size', 'mode'])['max'].agg(['median']).to_string())

# 191, 1s

# Plot modes for t1, c1 for each size
# for s in sorted(df['size'].unique()):
plt.figure()
df_t = df.copy()
df_t = df_t[df_t['threads'] == 1]
df_t = df_t[df_t['chunks'] == 1]
df_t = df_t[(df_t['size'] == 76800000) | (df_t['size'] == 191102976) | (df_t['size'] == 285153024)]
df_t = df_t.sort_values(['size'])
sns.set_style('darkgrid')
fig, axes = plt.subplots()
ax_violin = sns.violinplot(data = df_t, x = 'size', y = 'max', hue = 'mode', hue_order = ['datatype', 'manual'], ax = axes, facet_kws={'legend_out': True}, palette={'datatype':'tab:blue', 'manual':'tab:orange'}, inner=None)
a = sns.pointplot(data = df_t, x = 'size', y = 'max', hue = 'mode', hue_order = ['datatype', 'manual'], dodge = 0.4, ax = ax_violin, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.02)
plt.tick_params(axis='x', which='major', labelsize=9)
a.set_title("Redistribution time for manual packing and datatype", fontdict = {'fontsize' : 9})
plt.xlabel('Data size per node')
plt.ylabel("Time (μs)")
axes.set_xticklabels(['77MB (~13MB per send)', '191MB (~32MB per send)', '285MB (~48MB per send)'])
handles, labels = a.get_legend_handles_labels()
a.legend(handles, labels[:2], loc='upper left')
new_labels = ["Custom datatype", "Manual"]
for t, l in zip(a.get_legend().texts, new_labels):
    t.set_text(l)
plt.savefig(f'manual_vs_datatype.pdf')
plt.close()





# Plot threads for manual for every chunk, 191MB
plt.figure()
sns.set_style('darkgrid')
df_m = df.copy()
df_m = df_m[(df_m['mode'] == 'manual')]
df_m = df_m[df_m['chunks'] == 1]
df_m = df_m[df_m['size'] == 191102976]
df_m = df_m[df_m['threads'] <= 6]
fig, axes = plt.subplots()
cond = df_m['mode'] == 'datatype'
df_m.loc[cond,'threads'] = 0
ax_violin = sns.violinplot(data = df_m, x = 'threads', y = 'max', order=sorted(df_m['threads'].unique()), ax = axes, palette=['tab:orange', 'tab:orange', 'tab:orange', 'tab:orange', 'tab:orange', 'tab:orange'], inner=None)
a =  sns.pointplot(data = df_m, x = 'threads', y = 'max', order=sorted(df_m['threads'].unique()), ax = ax_violin, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.02)
plt.xlabel('Number of threads')
plt.ylabel("Time (μs)")
axes.set_xticklabels(['1', '2', '3', '4', '5', '6'])
a.set_title("Redistribution time (191MB per node) for manual packing with thread acceleration", fontdict = {'fontsize' : 9})
elems = a.get_children()
plt.legend(handles=[elems[0], elems[4]],labels=["Manual"], loc='upper left')
plt.savefig(f'manual_with_threads_191MB.pdf')
plt.close()




# Plot chunks for manual for every thread count, 191MB
plt.figure()
sns.set_style('darkgrid')
df_c = df.copy()
df_c = df_c[(df_c['mode'] == 'manual')]
df_c = df_c[df_c['threads'] == 1]
df_c = df_c[df_c['chunks'] <= 6]
df_c = df_c[df_c['size'] == 191102976]
fig, axes = plt.subplots()
cond = df_c['mode'] == 'datatype'
df_c.loc[cond,'chunks'] = 0
ax_violin = sns.violinplot(data = df_c, x = 'chunks', y = 'max', ax = axes, palette=['tab:orange', 'tab:orange', 'tab:orange', 'tab:orange', 'tab:orange', 'tab:orange'], inner=None)
a = sns.pointplot(data = df_c, x = 'chunks', y = 'max', order=sorted(df_c['chunks'].unique()), ax = ax_violin, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.02)
plt.xlabel('Number of chunks')
plt.ylabel("Time (μs)")
axes.set_xticklabels(['1', '2', '3', '4', '5', '6'])
a.set_title("Redistribution time (191MB per node) for manual packing with chunks", fontdict = {'fontsize' : 9})
elems = a.get_children()
plt.legend(handles=[elems[0], elems[4]],labels=["Manual"], loc='upper left')
plt.savefig(f'manual_packing_with_chunking_191MB.pdf')
plt.close()






















# # 191, 2s

# plt.figure()
# df_t = df.copy()
# df_t = df_t[df_t['threads'] == 1]
# fig, axes = plt.subplots()
# sns.violinplot(data = df_t, x = 'mode', y = 'max', ax = axes, inner=None)
# sns.pointplot(data = df_t, x = 'mode', y = 'max', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)
# plt.savefig(f'redist_modes_s191MB{file_suffix}.pdf')
# plt.close()

# plt.figure()
# df_m = df.copy()
# df_m = df_m[df_m['mode'] == 'manual']
# fig, axes = plt.subplots()
# sns.violinplot(data = df_m, x = 'threads', y = 'max', ax = axes, inner=None)
# sns.pointplot(data = df_m, x = 'threads', y = 'max', ax = axes, estimator=np.median, color='black', join=True, scale=0.5, errwidth=1, capsize=0.01)
# plt.savefig(f'redist_threads_s191MB{file_suffix}.pdf')
# plt.close()


# # Plot threads for manual for every chunk, 191MB plt.figure()
# df_m = df.copy()
# df_m = df_m[(df_m['mode'] == 'manual') | (df_m['mode'] == 'datatype')]
# df_m = df_m[df_m['chunks'] == 1]
# df_m = df_m[df_m['size'] == 191102976]
# fig, axes = plt.subplots()
# # sns.violinplot(data = df_m, x = 'threads', y = 'max', ax = axes, inner=None)
# sns.pointplot(data = df_m, x = 'threads', y = 'max', ax = axes, estimator=np.median, color='black', join=True, scale=0.5, errwidth=1, capsize=0.01)
# plt.savefig(f'redist_threads_s191102976_c1{file_suffix}.pdf')
# plt.close()

# # Plot chunks for manual for every thread count, 191MB
# plt.figure()
# df_c = df.copy()
# df_c = df_c[df_c['mode'] == 'manual']
# df_c = df_c[df_c['threads'] == 1]
# df_c = df_c[df_c['size'] == 191102976]
# fig, axes = plt.subplots()
# # sns.violinplot(data = df_c, x = 'chunks', y = 'max', ax = axes, inner=None)
# sns.pointplot(data = df_c, x = 'chunks', y = 'max', ax = axes, estimator=np.median, color='black', join=True, scale=0.5, errwidth=1, capsize=0.01)
# plt.savefig(f'redist_chunks_s191102976_t1{file_suffix}.pdf')
# plt.close()