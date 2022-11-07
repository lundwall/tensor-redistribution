import sys
import numpy as np

if len(sys.argv) != 2:
	print("exactly one input file should be specified.")
input_file = open(sys.argv[1], 'r')
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
		print(parameter_list)

statistics_array = np.array(statistics_list)
print(statistics_array)

