import matplotlib.pyplot as plt

threads = [1, 2, 3, 4, 5, 6, 7, 8]
medians = [31.054, 24.449, 32.292, 19.762, 18.264, 19.068, 21.310, 21.965]

ax = plt.gca()
ax.set_ylim([0, 35])

plt.plot(threads,medians)
plt.hlines(21.745, 1, 8, colors='red')
plt.text(1,20,'custom datatypes')
plt.title('5D Single Block Transmission Time by Number of Threads')
plt.xlabel('Threads #')
plt.ylabel('Time (ms)')

plt.savefig('multi_5d.pdf')

