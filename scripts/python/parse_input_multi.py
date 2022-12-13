import sys
import numpy as np
import os

if len(sys.argv) != 2:
	print("exactly one input folder should be specified.")

times = [[], [], [], [], [], [], [], []]

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
	print(statistics_array)

        

	times[int(statistics_array[0,1]) - 1] = times[int(statistics_array[0,1]) - 1] + statistics_array[:,3].tolist()
	print(times)


for t in range(8):
    print("With " + str(t + 1) + " threads")
    print("Mean: " + str(np.array([float(i) for i in times[t]]).mean()))
    print("STD: " + str(np.array([float(i) for i in times[t]]).std()))
    print("Median: " + str(np.median([float(i) for i in times[t]])))
    print("Min: " + str(np.array([float(i) for i in times[t]]).min()))
    print("Max: " + str(np.array([float(i) for i in times[t]]).max()))
    print("")

