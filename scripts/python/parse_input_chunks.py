from collections import defaultdict
import sys
import numpy as np
import os

if len(sys.argv) != 2:
	print("exactly one input folder should be specified.")

times = defaultdict(list)

for filename in os.listdir(os.getcwd() + '/' + sys.argv[1]):
	input_file = open(sys.argv[1] + '/' + filename, 'r')
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

        

	times[int(statistics_array[0,1])] += [float(e) for e in statistics_array[:,3].tolist()]
	# print(times)


with open('medians_chunks.csv', 'w') as f:
	f.write('Size,Time\n')
	for key, value in sorted(times.items()):
		value = np.array(value)
		print("With " + str(key) + " chunk size")
		print("Mean: " + str(value.mean()))
		print("STD: " + str(value.std()))
		print("Median: " + str(np.median(value)))
		print("Min: " + str(value.min()))
		print("Max: " + str(value.max()))
		print("")
		f.write(str(key) + ',' + str(np.median(value)) + '\n')


