import numpy as np
import matplotlib.pyplot as plt
from mpl_point_clicker import clicker

fig, ax = plt.subplots()
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
klicker = clicker(ax, ["event"], markers=["o"], **{"linestyle": "--"})

plt.show()

print('[')
for coord in klicker.get_positions()['event']:
    print(' [', coord[0], ', ', coord[1], '],\n', sep='', end='')
print(']', end='')