import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imageio.v2 as imageio
import os
import random

# Ensure the frames directory exists
os.makedirs('./frames', exist_ok=True)

# Open the HDF5 file and read the dataset
with h5py.File('electron_positions.h5', 'r', encoding='utf-8') as f:
    frames = f['positions'][:]  # Read all the frames into memory
# List of 10 color options
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown']

# Assign a random color to each point in the first frame
point_colors = [random.choice(colors) for _ in range(frames.shape[1])]

# Initialize previous frame as None
previous_frame = None
image_files = []
# Plot each frame
for i in range(0, len(frames), 10):
    frame = frames[i]
    fig = plt.figure(figsize=(10, 8))  # Adjust the figure size
    ax = fig.add_subplot(111, projection='3d')
    x, y, z = frame[:, 0], frame[:, 1], frame[:, 2]

    # Use the assigned colors for the points
    ax.scatter(x, y, z, color=point_colors)
    # If there is a previous frame, draw lines between the points
    if previous_frame is not None:
        prev_x, prev_y, prev_z = previous_frame[:, 0], previous_frame[:, 1], previous_frame[:, 2]
        for j in range(len(frame)):
            ax.plot([prev_x[j], x[j]], [prev_y[j], y[j]], [prev_z[j], z[j]], color='r')
    ax.set_title(f'Time={i * 0.0001:.4f}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.set_zlim(-10, 10)
    plt.grid(True)
    filename = f'./frames/frame_{i}.png'
    plt.savefig(f'./frames/frame_{i}.png')
    plt.close()
    image_files.append(filename)
    previous_frame = frame

# List to store the file names of the generated images


# Create a list to store the images
images = []

# Save the images as a GIF
with imageio.get_writer('animation.gif', mode='I') as writer:
    for filename in image_files:
        img = imageio.imread(filename)
        writer.append_data(img, {'duration': 0.1})