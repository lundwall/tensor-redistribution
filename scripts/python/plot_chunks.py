import matplotlib.pyplot as plt
import pandas as pd

columns = ["Size", "Time"]
df = pd.read_csv("medians_chunks.csv", usecols=columns)

ax = plt.gca()
ax.set_xscale('log')

# start_index = 18
# plt.plot(df.Size.to_list()[start_index:], df.Time.to_list()[start_index:])
plt.plot(df.Size, df.Time)
plt.title('2D Single Block Transmission Time by Chunk Size')
plt.xlabel('Chunk Size (ints)')
plt.ylabel('Median Time (ms)')

plt.savefig('chunk_2d_log_split.pdf')

