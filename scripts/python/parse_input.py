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

statistics_array = np.array(statistics_list)

times = statistics_array[:,3]

desc = ''
if '2d' in sys.argv[1]:
	desc = desc + '2D '
elif '3d' in sys.argv[1]:
	desc = desc + '3D '
elif '4d' in sys.argv[1]:
	desc = desc + '4D '
if 'manual' in sys.argv[1]:
	desc = desc + 'with manual packing, '
elif 'datatype' in sys.argv[1]:
	desc = desc + 'with custom datatypes, '
if 'r0' in sys.argv[1]:
	desc = desc + 'rank 0'
elif 'r1' in sys.argv[1]:
	desc = desc + 'rank 1'

print(desc)
print("Mean: " + str(np.array([float(i) for i in times]).mean()))
print("STD: " + str(np.array([float(i) for i in times]).std()))
print("Median: " + str(np.median([float(i) for i in times])))
print("Min: " + str(np.array([float(i) for i in times]).min()))
print("Max: " + str(np.array([float(i) for i in times]).max()))
print("")

