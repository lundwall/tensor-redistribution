import sys
import numpy as np
import os

if len(sys.argv) != 2:
	print("exactly one input folder should be specified.")

times = []

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

	times = times + statistics_array[:,4].tolist()

print("Mean: " + str(np.array([float(i) for i in times]).mean()))
print("STD: " + str(np.array([float(i) for i in times]).std()))
print("Median: " + str(np.median([float(i) for i in times])))
print("Min: " + str(np.array([float(i) for i in times]).min()))
print("Max: " + str(np.array([float(i) for i in times]).max()))
print("")

