import numpy as np
import matplotlib.pyplot as plt


#load output.txt as list of strings
with open('output.txt') as f:
    lines = f.readlines()
#remove newline characters
lines = [line.strip() for line in lines]

#separate by "___"


fieldStrings = [[]]
for line in lines:
    if "_" in line:
        fieldStrings.append([])
    else:
        fieldStrings[-1].append(line)

#print number of field Strings
print("NT:")
print(len(fieldStrings))

#now parse into numbers (each line is csv)
fieldData = []
for fieldString in fieldStrings:
    timeStepData = []
    for line in fieldString:
        timeStepData.append([float(x) for x in line.split(",")])
    fieldData.append(timeStepData)

#print number of floats in the first row of the first field
print(len(fieldData[0][0]))

#plot the first field as a 2D image
fieldData1 = np.array(fieldData[99])
plt.imshow(fieldData1)

#save plot to .png
plt.savefig('fieldData1.png')

import matplotlib.animation as animation

# Assuming fieldData is already parsed and loaded as a 3D list (NT, NX, NY)

# Set up the figure, the axis, and the plot element we want to animate
fig, ax = plt.subplots()

# Initialize the plot with the first frame
cax = ax.matshow(fieldData[0], cmap='viridis')
fig.colorbar(cax)

def update(frame):
    ax.clear()  # Clear the previous frame
    cax = ax.matshow(fieldData[frame], cmap='viridis')  # Plot the new frame
    ax.set_title(f'Time Step {frame+1}')
    return cax,

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=len(fieldData), blit=False)

# Save the animation as a GIF using Pillow writer
ani.save("output_animation.gif", writer='pillow', fps=20)

print("Animation saved as output_animation.gif")