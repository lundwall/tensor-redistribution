import sys
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
	print("exactly one input folder should be specified.")

times = [[], [], [], [], [], [], [], []]
cases = ['5d']
node_types = ['r0']
dict = [{}, {}, {}, {}, {}, {}, {}, {}]

experiments = ['with_API', 'without_API', 'manual_put', 'custom_datatype', 'custom_datatype_put']
for rank in node_types:
	for dim in cases:
		for num_thread in range(1, 9):
			for experiment in experiments:
				filename = 'lsb.' + dim + '_transmit_' + experiment +str(num_thread)+'.' + rank
				input_file = open(sys.argv[1]+'/' + filename, 'r')
				lines = input_file.readlines()
				parameter_list = []
				statistics_list = []

				for line in lines:
					elements = line.split()
					if line[0] == '#':
						continue
					try:
						_ = float(elements[0])
						statistics_list.append(elements)
					except ValueError:
						parameter_list = elements

				statistics_array = np.array(statistics_list)
				if experiment == 'with_API':
					datatype_time = statistics_array[:, 7]
				else:
					datatype_time = statistics_array[:, 2]
				median = np.median([float(i) for i in datatype_time])
				max_num = np.max([float(i) for i in datatype_time])
				new_time = []
				for t in datatype_time:
					if float(t) < median * 15:
						new_time.append(t)
				dict[num_thread-1][filename] = new_time

sns.set(font_scale=0.6)

# for i in range(8):
# 	fig, ax = plt.subplots(5, 1)
# 	fig.set_size_inches(6, 8)
# 	fig.subplots_adjust(hspace=0.5)
# 	for key, ax_i in zip(dict[i].keys(), ax.reshape(-1)):
# 		sns.distplot(dict[i][key], norm_hist=True, bins=30, ax=ax_i).set(title=key)
# 	plt.savefig('different_method_same_threadnum_' + str(i+1) + '.pdf')


fig, ax = plt.subplots(3, 3)
fig.set_size_inches(8, 8)
fig.subplots_adjust(hspace=0.5, wspace=0.5, left=0.1, right=0.9)
for num_thread, ax_i in zip(range(1, 9), ax.reshape(-1)):
	filename = 'lsb.5d_transmit_without_API' + str(num_thread) +'.r0'
	sns.distplot(dict[num_thread-1][filename], bins=30, ax=ax_i).set(title=filename)
	filename = 'lsb.5d_transmit_custom_datatype1.r0'
	sns.distplot(dict[0][filename], bins=30, ax=ax_i)

filename = 'lsb.5d_transmit_custom_datatype1.r0'
sns.distplot(dict[0][filename],bins=30, ax=ax[2,2]).set(title='custom_datatype')
plt.savefig('manual_vs_datatype.pdf')



			# for num_thread in range(1, 9):
			# # for filename in os.listdir(os.getcwd()):
			# 	filename = 'lsb.' + dim + '_transmit_manual' + str(num_thread) + '.' + rank
			# 	# print(filename)
			# 	input_file = open(sys.argv[1]+'/' + filename, 'r')
			# 	lines = input_file.readlines()
			# 	parameter_list = []
			# 	statistics_list = []
			#
			# 	for line in lines:
			# 		elements = line.split()
			# 		if line[0] == '#':
			# 			continue
			# 		try:
			# 			_ = float(elements[0])
			# 			statistics_list.append(elements)
			# 		except ValueError:
			# 			parameter_list = elements
			#
			# 	statistics_array = np.array(statistics_list)
			# 	# print(statistics_array)
			#
			# 	times[int(statistics_array[0, 1]) - 1] = times[int(statistics_array[0, 1]) - 1] + statistics_array[:,
			# 																					  3].tolist()
			# # print(times)
			#
			# Medians = []
			# Means = []
			# for t in range(8):
			# 	# print("With " + str(t + 1) + " threads")
			# 	# print("Mean: " + str(np.array([float(i) for i in times[t]]).mean()))
			# 	Means.append(np.array([float(i) for i in times[t]]).mean() / 1000)
			# 	# print("STD: " + str(np.array([float(i) for i in times[t]]).std()))
			# 	# print("Median: " + str(np.median([float(i) for i in times[t]])))
			# 	Medians.append(np.median([float(i) for i in times[t]]) / 1000)
			# # print("Min: " + str(np.array([float(i) for i in times[t]]).min()))
			# # print("Max: " + str(np.array([float(i) for i in times[t]]).max()))
			# # print("")
			#
			# # print(dim, "manual means:", Means)
			# # print("\n")
			# print(rank, dim, "manual medians:", Medians)
			# print("\n")

