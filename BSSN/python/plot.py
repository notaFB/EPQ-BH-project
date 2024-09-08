# Load necessary libraries
import numpy as np
import matplotlib.pyplot as plt

# File paths
path_phi = "/home/andres/Desktop/not/NOT_A/EPQ-BH-project/BSSN/build/phi_data.txt"
path_phi_pre = "/home/andres/Desktop/not/NOT_A/EPQ-BH-project/BSSN/build/phi_pre_data.txt"

# Load the data from the files
data_phi = np.loadtxt(path_phi)
data_phi_pre = np.loadtxt(path_phi_pre)

# Reshape the data into a 3D grid, since N=81, we reshape into (N, N, N)
N = 31
data_phi_3d = data_phi.reshape((N, N, N))
data_phi_pre_3d = data_phi_pre.reshape((N, N, N))

# Select the middle slice for z (z index = 40)
z_middle = N // 2
middle_slice_phi = data_phi_3d[:, :, z_middle]
middle_slice_phi_pre = data_phi_pre_3d[:, :, z_middle]

# Create subplots to display both datasets side by side
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# Plot the first dataset (phi_data)
im1 = axes[0].imshow(middle_slice_phi, cmap='viridis', origin='lower')
axes[0].set_title(f"Phi Data - Middle Z Slice (z index = {z_middle})")
axes[0].set_xlabel("X Index")
axes[0].set_ylabel("Y Index")
fig.colorbar(im1, ax=axes[0], label="Phi Value")

# Plot the second dataset (phi_pre_data)
im2 = axes[1].imshow(middle_slice_phi_pre, cmap='viridis', origin='lower')
axes[1].set_title(f"Phi Pre Data - Middle Z Slice (z index = {z_middle})")
axes[1].set_xlabel("X Index")
axes[1].set_ylabel("Y Index")
fig.colorbar(im2, ax=axes[1], label="Phi Pre Value")

# Show the plots
plt.tight_layout()
plt.show()
