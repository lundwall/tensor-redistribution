from matplotlib import ticker as mticker
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
        # sometimes like (12,12,12,12,96) has better chunking result then (24,24,24,24,24). So when running
        # experiments I also include block size like that. But I'm not sure whether they should be put inside 
        # the result plots. Comment out the following if statement if you want to display the size also.
        if size_tuple[0] != size_tuple[1] or size_tuple[1] != size_tuple[2] or size_tuple[2] != size_tuple[3] or size_tuple[3] != size_tuple[4]:
            continue
        if size_tuple in size_to_df:
            size_to_df[size_tuple].append(df)
        else:
            size_to_df[size_tuple] = [df]
    except ValueError:
        print("some size is not integer")

total_df = []
for key, value in size_to_df.items():
    size_to_df[key] = pd.concat(value)
    total_size = key[0] * key[1] * key[2] * key[3] * key[4] * 4 / 1000000
    # 12**5 is too small, skip it in processing
    if total_size < 1:
        continue
    append_column = [total_size for x in range(len(size_to_df[key].index))]
    size_to_df[key] = size_to_df[key].assign(total_size=append_column)

    size_string = str(key[0]) + "*" + str(key[1]) + "*" + str(key[2]) + "*" + str(key[3]) + "*" + str(key[4])
    append_column = [size_string for x in range(len(size_to_df[key].index))]
    size_to_df[key] = size_to_df[key].assign(size_string=append_column)
    
    total_df.append(size_to_df[key])

df_all = pd.concat(total_df)
df_all = df_all[df_all['type'] == 'running'] # remove syncing runs

# manual packing
df_t1_c1 = df_all[(df_all['chunk'] == 1) & (df_all['threads'] == 1)]
df_t1_c1 = df_t1_c1[(df_t1_c1['mode'] == 'datatype') | (df_t1_c1['mode'] == 'manual')]
df_t1_c1 = df_t1_c1.groupby(['id', 'mode', 'size_string', 'total_size'])['time'].agg(['max']).reset_index()
df_t1_c1 = df_t1_c1.sort_values(['total_size'])

plt.figure()
sns.set_style('darkgrid')
fig, axes = plt.subplots()
df_t1_c1['max'] =  np.log2(df_t1_c1['max'])
a = sns.violinplot(data = df_t1_c1, x = 'size_string', y = 'max', hue = 'mode', hue_order = ['datatype', 'manual'], ax = axes, facet_kws={'legend_out': True})
# sns.pointplot(data = df_t1_c1, x = 'mode', y = 'max', oder=sorted(df['mode'].unique()), ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)
plt.tick_params(axis='x', which='major', labelsize=9)
a.set_title("Single-block transmission time for manual packing and datatype", fontdict = {'fontsize' : 9})
plt.xlabel('Block size')
plt.ylabel("Time(us)")

axes.yaxis.set_major_formatter(mticker.StrMethodFormatter("$2^{{{x:.0f}}}$"))
ymin, ymax = axes.get_ylim()
tick_range = np.arange(np.floor(ymin), ymax)
size_list = sorted(df_t1_c1['total_size'].unique())
axes.set_xticklabels([str(int(x))+'MB' for x in size_list])
axes.yaxis.set_ticks(tick_range)
new_labels = ['datatype', 'manual']
for t, l in zip(axes.get_legend().texts, new_labels):
    t.set_text(l)
plt.savefig(f'manual_vs_datatype.pdf')
plt.close()

# manual packing with thread acceleration 47MB
df_c1_lsize = df_all[((df_all['size_string'] == '26*26*26*26*26') & (df_all['mode'] == 'manual') & (df_all['chunk'] == 1)) | ((df_all['size_string'] == '26*26*26*26*26') & (df_all['mode'] == 'datatype') & (df_all['threads'] == 1))]
df_c1_lsize = df_c1_lsize.groupby(['id', 'mode', 'threads'])['time'].agg(['max']).reset_index()
df_c1_lsize['max'] =  np.log2(df_c1_lsize['max'])
plt.figure()
sns.set_style('darkgrid')
fig, axes = plt.subplots()
# hack: set datatype thread = 0 here so it would stay at left most part of the plot
cond = df_c1_lsize['mode'] == 'datatype'
df_c1_lsize.loc[cond,'threads'] = 0
a = sns.violinplot(data = df_c1_lsize, x = 'threads', y = 'max', order=sorted(df_c1_lsize['threads'].unique()),ax = axes)
# sns.pointplot(data = df_c1_lsize, x = 'chunk', y = 'max', order=sorted(df_c1_lsize['threads'].unique()), hue = 'mode', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)
plt.xlabel('Thread number')
plt.ylabel('Time(us)')
axes.set_xticklabels(['datatype', '1', '2', '3', '4'])
axes.yaxis.set_major_formatter(mticker.StrMethodFormatter("$2^{{{x:.0f}}}$"))
ymin, ymax = axes.get_ylim()
tick_range = np.arange(np.floor(ymin), ymax)
axes.yaxis.set_ticks(tick_range)
a.set_title("Single-block(47MB) transmission time for manual packing with thread acceleration", fontdict = {'fontsize' : 9})
plt.savefig(f'manual_with_threads_vs_datatype_47MB.pdf')
plt.close()

# manual packing with thread acceleration 31MB
df_c1_lsize = df_all[((df_all['size_string'] == '24*24*24*24*24') & (df_all['mode'] == 'manual') & (df_all['chunk'] == 1)) | ((df_all['size_string'] == '24*24*24*24*24') & (df_all['mode'] == 'datatype') & (df_all['threads'] == 1))]
df_c1_lsize = df_c1_lsize.groupby(['id', 'mode', 'threads'])['time'].agg(['max']).reset_index()
df_c1_lsize['max'] =  np.log2(df_c1_lsize['max'])
plt.figure()
sns.set_style('darkgrid')
fig, axes = plt.subplots()
# hack: set datatype thread = 0 here so it would stay at left most part of the plot
cond = df_c1_lsize['mode'] == 'datatype'
df_c1_lsize.loc[cond,'threads'] = 0
a = sns.violinplot(data = df_c1_lsize, x = 'threads', y = 'max', order=sorted(df_c1_lsize['threads'].unique()), ax = axes)
# sns.pointplot(data = df_c1_lsize, x = 'chunk', y = 'max', order=sorted(df_c1_lsize['threads'].unique()), hue = 'mode', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)
plt.xlabel('Thread number')
plt.ylabel("Time(us)")
axes.set_xticklabels(['datatype', '1', '2', '3', '4'])
axes.yaxis.set_major_formatter(mticker.StrMethodFormatter("$2^{{{x:.0f}}}$"))
ymin, ymax = axes.get_ylim()
tick_range = np.arange(np.floor(ymin), ymax)
axes.yaxis.set_ticks(tick_range)
a.set_title("Single-block(31MB) transmission time for manual packing with thread acceleration", fontdict = {'fontsize' : 9})
plt.savefig(f'manual_with_threads_vs_datatype_31MB.pdf')
plt.close()

# manual packing with chunking 47MB
df_t1_lsize = df_all[(df_all['size_string'] == '26*26*26*26*26') & (df_all['threads'] == 1)]
df_t1_lsize = df_t1_lsize.groupby(['id', 'mode', 'chunk'])['time'].agg(['max']).reset_index()
df_t1_lsize = df_t1_lsize[(df_t1_lsize['mode'] == 'datatype') | (df_t1_lsize['mode'] == 'manual')]
df_t1_lsize['max'] =  np.log2(df_t1_lsize['max'])
plt.figure()
# hack: set datatype chunk = 0 here so it would stay at left most part of the plot
cond = df_t1_lsize['mode'] == 'datatype'
df_t1_lsize.loc[cond,'chunk'] = 0
sns.set_style('darkgrid')
fig, axes = plt.subplots()
a = sns.violinplot(data = df_t1_lsize, x = 'chunk', y = 'max', ax = axes)
# sns.pointplot(data = df_t1_lsize, x = 'chunk', y = 'max', order=sorted(df_t1_lsize['chunk'].unique()), hue = 'mode', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)
plt.xlabel('Chunk number')
plt.ylabel("Time(us)")
axes.set_xticklabels(['datatype', '1', '2', '3', '4', '5'])
axes.yaxis.set_major_formatter(mticker.StrMethodFormatter("$2^{{{x:.0f}}}$"))
ymin, ymax = axes.get_ylim()
tick_range = np.arange(np.floor(ymin), ymax)
axes.yaxis.set_ticks(tick_range)
a.set_title("Single-block(47MB) transmission time for manual packing with chunks", fontdict = {'fontsize' : 9})
plt.savefig(f'manual_packing_with_chunking_vs_datatype_47MB.pdf')
plt.close()

# manual packing with chunking 31MB
df_t1_lsize = df_all[(df_all['size_string'] == '24*24*24*24*24') & (df_all['threads'] == 1)]
df_t1_lsize = df_t1_lsize.groupby(['id', 'mode', 'chunk'])['time'].agg(['max']).reset_index()
df_t1_lsize = df_t1_lsize[(df_t1_lsize['mode'] == 'datatype') | (df_t1_lsize['mode'] == 'manual')]
df_t1_lsize['max'] =  np.log2(df_t1_lsize['max'])
plt.figure()
# hack: set datatype chunk = 0 here so it would stay at left most part of the plot
cond = df_t1_lsize['mode'] == 'datatype'
df_t1_lsize.loc[cond,'chunk'] = 0
sns.set_style('darkgrid')
fig, axes = plt.subplots()
a = sns.violinplot(data = df_t1_lsize, x = 'chunk', y = 'max', ax = axes)
# sns.pointplot(data = df_t1_lsize, x = 'chunk', y = 'max', order=sorted(df_t1_lsize['chunk'].unique()), hue = 'mode', ax = axes, estimator=np.median, color='black', join=False, scale=0.5, errwidth=1, capsize=0.01)
plt.xlabel('Chunk number')
plt.ylabel("Time(us)")
axes.set_xticklabels(['datatype', '1', '2', '3', '4', '5'])
axes.yaxis.set_major_formatter(mticker.StrMethodFormatter("$2^{{{x:.0f}}}$"))
ymin, ymax = axes.get_ylim()
tick_range = np.arange(np.floor(ymin), ymax)
axes.yaxis.set_ticks(tick_range)
a.set_title("Single-block(31MB) transmission time for manual packing with chunks", fontdict = {'fontsize' : 9})
plt.savefig(f'manual_packing_with_chunking_vs_datatype_31MB.pdf')
plt.close()

# MPI one sided
df_c1_lsize = df_all[(df_all['chunk'] == 1) & (df_all['threads'] == 1)]
df_c1_lsize = df_c1_lsize.groupby(['id', 'mode', 'size_string', 'total_size'])['time'].agg(['max']).reset_index()
df_c1_lsize = df_c1_lsize.sort_values(['total_size'])
df_c1_lsize['max'] =  np.log2(df_c1_lsize['max'])
plt.figure()
sns.set_style('darkgrid')
fig, axes = plt.subplots()
a = sns.violinplot(data = df_c1_lsize, x = 'size_string', y = 'max', hue = 'mode', hue_order = ['datatype', 'manual', 'put_manual'], ax = axes, facet_kws={'legend_out': True})
plt.tick_params(axis='x', which='major', labelsize=9)
plt.xlabel('Block size')
plt.ylabel("Time(us)")
size_list = sorted(df_c1_lsize['total_size'].unique())
axes.set_xticklabels([str(int(x))+'MB' for x in size_list])
axes.yaxis.set_major_formatter(mticker.StrMethodFormatter("$2^{{{x:.0f}}}$"))
ymin, ymax = axes.get_ylim()
tick_range = np.arange(np.floor(ymin), ymax)
axes.yaxis.set_ticks(tick_range)
a.set_title("Single-block transmission time for datatype, one-sided and manual packing", fontdict = {'fontsize' : 9})
new_labels = ['datatype', 'manual', 'one-sided']
for t, l in zip(axes.get_legend().texts, new_labels):
    t.set_text(l)
plt.savefig(f'manual_packing_with_one_sided_vs_datatype.pdf')
plt.close()