import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import matplotlib.animation as animation
# from matplotlib.animation import PillowWriter

N = 200
K = 1.0
rho = 0.2
L = np.sqrt(N/rho)

data_folder = "movie_data/"
file_list = sorted(glob.glob(os.path.join(data_folder, "step_*.txt")))

# file_list = file_list[-1000:]

data = np.loadtxt(file_list[0])
x, y, theta = data[:, 0], data[:, 1], data[:, 2]


color_values = np.cos(theta)
fig, ax = plt.subplots()
sc = ax.scatter(x, y, c=color_values, cmap='RdBu', vmin=-1, vmax=1, s=10)
ax.set_xlim(0, L+1)
ax.set_ylim(0, L+1)


cbar = plt.colorbar(sc)
cbar.set_label(r'$\cos(\theta_i)$')

plt.gca().set_box_aspect(1)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("AH Model Evolution")


def update(frame):
    data = np.loadtxt(file_list[frame])
    x, y, theta = data[:, 0], data[:, 1], data[:, 2]
    color_values = np.cos(theta)
    sc.set_offsets(np.c_[x, y])  # Update particle positions
    sc.set_array(color_values)    # Update colors based on theta
    ax.set_title(f"AH Model - Step {frame}")

# Create animation
ani = animation.FuncAnimation(fig, update, frames=len(file_list), interval=5, repeat=True)
# ani.save("ah_system_anim.gif", writer=PillowWriter(fps=60))
# Show animation
plt.show()

