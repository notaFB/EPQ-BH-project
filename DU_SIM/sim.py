import numpy as np
import sympy as sp
import mpmath
import pandas as pd
import matplotlib.pyplot as plt

# Set precision for mpmath
mpmath.mp.dps = 50  # Set to 50 decimal places for higher precision

# 1. Define the Minkowski metric in spherical coordinates
t, r, theta, phi = sp.symbols('t r theta phi')
coords_spherical = [t, r, theta, phi]

g_spherical = sp.Matrix([
    [-1, 0, 0, 0],
    [ 0, 1, 0, 0],
    [ 0, 0, r**2, 0],
    [ 0, 0, 0, r**2 * sp.sin(theta)**2]
])

g_spherical_inv = g_spherical.inv()

# Define the Christoffel symbol function
def GammaFromMetric(g, ginv, coords):
    dg = sp.derive_by_array(g, coords)
    g_ab_d = sp.permutedims(dg.copy(), (1, 2, 0))
    g_bd_a = sp.permutedims(dg.copy(), (0, 1, 2))
    g_ad_b = sp.permutedims(dg.copy(), (1, 0, 2))
    partials_g = g_ad_b + g_bd_a - g_ab_d
    christoffel = sp.tensorcontraction(sp.tensorproduct(partials_g, ginv), (2, 4)) * 0.5
    christoffel = sp.permutedims(christoffel, (2, 0, 1))
    return christoffel

christoffel_spherical = GammaFromMetric(g_spherical, g_spherical_inv, coords_spherical)

# 2. Lambdify the Christoffel symbols for numerical evaluation using mpmath for higher precision
christoffel_func = sp.lambdify((t, r, theta, phi), christoffel_spherical, 'mpmath')

# 3. Define the geodesic equation using correct indices and loops
def geodesic_equation(X, christoffel_func):
    # X contains (t, r, theta, phi, dtdtau, drdtau, dthetadtau, dphidtau)
    t_val, r_val, theta_val, phi_val = X[0], X[1], X[2], X[3]
    derivatives = X[4:]  # (dtdtau, drdtau, dthetadtau, dphidtau)
    
    Gamma = np.array(christoffel_func(t_val, r_val, theta_val, phi_val), dtype=mpmath.mpf)

    # Initialize second derivatives (d2x/dtau^2 for each coordinate)
    second_derivatives = np.zeros(4, dtype=mpmath.mpf)
    
    for mu in range(4):  # Loop over each coordinate (t, r, theta, phi)
        for alpha in range(4):
            for beta in range(4):
                second_derivatives[mu] += -Gamma[mu, alpha, beta] * derivatives[alpha] * derivatives[beta]
    
    return np.concatenate([derivatives, second_derivatives])  # Returns [dx/dtau, d2x/dtau^2]

# RK4 integrator using mpmath
def rk4_step(X, christoffel_func, dtau):
    k1 = geodesic_equation(X, christoffel_func)
    k2 = geodesic_equation(X + 0.5 * dtau * k1, christoffel_func)
    k3 = geodesic_equation(X + 0.5 * dtau * k2, christoffel_func)
    k4 = geodesic_equation(X + dtau * k3, christoffel_func)
    return X + (dtau / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

# 3.1 Initial conditions for a straight line (with mpmath precision)
#X0 -> [t, r, theta, phi, dtdtau, drdtau, dthetadtau, dphidtau]
X0_straight = np.array([mpmath.mpf(3.0), mpmath.mpf(1.0), mpmath.pi / 2, mpmath.mpf(0.0), mpmath.mpf(1.0), mpmath.mpf(0.3), mpmath.mpf(0.0), mpmath.mpf(1.0)])

# Integration parameters
tau_start = mpmath.mpf(0.0)
tau_end = mpmath.mpf(2.0)
dtau = mpmath.mpf(0.01)  # Smaller step size for higher precision
N_steps = int((tau_end - tau_start) / dtau)

# Define precision arrays to store results as we go
tau_vals = []
t_vals = []
r_vals = []
theta_vals = []
phi_vals = []

# Modify the RK4 solver loop to store values in memory during iterations
X = X0_straight
for i in range(N_steps):
    tau_vals.append(float(tau_start + i * dtau))  # Store tau value
    t_vals.append(X[0])
    r_vals.append(X[1])
    theta_vals.append(mpmath.pi / 2)  # Constant theta
    phi_vals.append(X[3])
    
    # Take the next RK4 step
    X = rk4_step(X, christoffel_func, dtau)

# Now let's save the data into a CSV file after the loop
trajectory_data = {
    'tau': tau_vals,
    't': t_vals,
    'r': r_vals,
    'theta': theta_vals,
    'phi': phi_vals
}

# Convert the data into a pandas DataFrame for easier handling
trajectory_df = pd.DataFrame(trajectory_data)

# Save the DataFrame to a CSV file with full precision
file_path = './simdata.csv'
trajectory_df.to_csv(file_path, float_format='%.50f', index=False)

# 4. Convert from spherical (r, phi) -> Cartesian (x, y) using mpmath precision
x_vals = np.array([r * mpmath.cos(phi) for r, phi in zip(r_vals, phi_vals)], dtype=mpmath.mpf)
y_vals = np.array([r * mpmath.sin(phi) for r, phi in zip(r_vals, phi_vals)], dtype=mpmath.mpf)

# 4. Plot the trajectory in 2D using float conversion for plotting
plt.figure(figsize=(8, 6))
plt.plot([float(x) for x in x_vals], [float(y) for y in y_vals], label='Geodesic trajectory (expected straight line)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('2D Geodesic Plot (r, phi) -> (x, y) with Correct Christoffel Symbols (dtau = 0.01)')
plt.legend()
plt.grid(True)
#ranges for x 0, 4 and y 0, 2
plt.xlim(0, 4)
plt.ylim(0, 2)
plt.show()


