import sys
import numpy as np
import os


times = [[], [], [], [], [], [], [], []]
cases = ['2d', '3d', '4d', '5d']
node_types = ['r0', 'r1']
for rank in node_types:
	for dim in cases:
		filename = 'lsb.' + dim + '_transmit_datatype.' + rank
		input_file = open('./' + filename, 'r')
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
		# print(statistics_array)
		datatype_time = statistics_array[:,3]
		print(rank, dim, 'datatype median', np.median([float(i) for i in datatype_time]) / 1000)

		for num_thread in range(1, 9):
		# for filename in os.listdir(os.getcwd()):
			filename = 'lsb.' + dim + '_transmit_manual' + str(num_thread) + '.' + rank
			# print(filename)
			input_file = open('./' + filename, 'r')
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
			# print(statistics_array)

			times[int(statistics_array[0, 1]) - 1] = times[int(statistics_array[0, 1]) - 1] + statistics_array[:,
																							  3].tolist()
		# print(times)

		Medians = []
		Means = []
		for t in range(8):
			# print("With " + str(t + 1) + " threads")
			# print("Mean: " + str(np.array([float(i) for i in times[t]]).mean()))
			Means.append(np.array([float(i) for i in times[t]]).mean() / 1000)
			# print("STD: " + str(np.array([float(i) for i in times[t]]).std()))
			# print("Median: " + str(np.median([float(i) for i in times[t]])))
			Medians.append(np.median([float(i) for i in times[t]]) / 1000)
		# print("Min: " + str(np.array([float(i) for i in times[t]]).min()))
		# print("Max: " + str(np.array([float(i) for i in times[t]]).max()))
		# print("")

		# print(dim, "manual means:", Means)
		# print("\n")
		print(rank, dim, "manual medians:", Medians)
		print("\n")

