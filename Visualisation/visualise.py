import json
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imageio.v2 as imageio
import os
import random

# Ensure the frames directory exists
os.makedirs('./frames', exist_ok=True)

# Read frames from JSON file
with open('frames.json', 'r') as f:
    frames = json.load(f)

# List of 10 color options
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown']

# Assign a random color to each point in the first frame
point_colors = [random.choice(colors) for _ in frames[0]]

# Initialize previous frame as None
previous_frame = None

# Plot each frame
for i, frame in enumerate(frames):
    fig = plt.figure(figsize=(10, 8))  # Adjust the figure size
    ax = fig.add_subplot(111, projection='3d')
    x, y, z = zip(*frame)
    
    # Use the assigned colors for the points
    ax.scatter(x, y, z, color=point_colors)
    
    # If there is a previous frame, draw lines between the points
    if previous_frame is not None:
        prev_x, prev_y, prev_z = zip(*previous_frame)
        for j in range(len(frame)):
            ax.plot([prev_x[j], x[j]], [prev_y[j], y[j]], [prev_z[j], z[j]], color='r')
    
    ax.set_title(f'Time={i * 0.01:.2f}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.set_zlim(-10, 10)
    plt.grid(True)
    plt.savefig(f'./frames/frame_{i}.png')
    plt.close()
    
    # Update previous frame
    previous_frame = frame

# List to store the file names of the generated images
image_files = [f'./frames/frame_{i}.png' for i in range(len(frames))]

# Read images and store them in a list
images = []
for filename in image_files:
    images.append(imageio.imread(filename))

# Save the images as a GIF
imageio.mimsave('animation.gif', images, duration=0.5)