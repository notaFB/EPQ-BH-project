#include <cmath>
#include <iostream>

// Grid dimensions
const int NX = 100; // Number of grid points in x-direction
const int NY = 100; // Number of grid points in y-direction
const int NZ = 100; // Number of grid points in z-direction

// Physical grid spacing
const double dx = 1.0 / NX;
const double dy = 1.0 / NY;
const double dz = 1.0 / NZ;

// BSSN variables
double phi[NX][NY][NZ];                 // Conformal factor
double K[NX][NY][NZ];                   // Trace of the extrinsic curvature
double gamma_tilde[NX][NY][NZ][3][3];   // Conformal metric
double A_tilde[NX][NY][NZ][3][3];       // Traceless part of the conformal extrinsic curvature
double Gamma_tilde[NX][NY][NZ][3];      // Conformal connection functions

// Gauge variables
double alpha[NX][NY][NZ];               // Lapse function
double beta[NX][NY][NZ][3];             // Shift vector
double B[NX][NY][NZ][3];                // Auxiliary variable (for Gamma driver)

void initialize_grid() {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                // Set the conformal metric to be flat
                gamma_tilde[i][j][k][0][0] = 1.0;
                gamma_tilde[i][j][k][1][1] = 1.0;
                gamma_tilde[i][j][k][2][2] = 1.0;
                gamma_tilde[i][j][k][0][1] = 0.0;
                gamma_tilde[i][j][k][0][2] = 0.0;
                gamma_tilde[i][j][k][1][2] = 0.0;
                
                // Initialize the conformal factor
                phi[i][j][k] = 1.0; // Start with a uniform guess

                // Initialize the extrinsic curvature components to zero
                K[i][j][k] = 0.0;

                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        A_tilde[i][j][k][l][m] = 0.0;
                    }
                    Gamma_tilde[i][j][k][l] = 0.0;
                }

                // Initialize lapse and shift
                alpha[i][j][k] = 1.0;  // Start with alpha = 1 everywhere
                for (int l = 0; l < 3; l++) {
                    beta[i][j][k][l] = 0.0;  // Start with beta = 0
                    B[i][j][k][l] = 0.0;     // Initialize B to zero
                }
            }
        }
    }
}

// Constants for the black holes
const double M1 = 0.5;  // Mass of the first black hole
const double M2 = 0.5;  // Mass of the second black hole
const double P1 = 0.1;  // Momentum of the first black hole
const double P2 = -0.1; // Momentum of the second black hole
const double X1 = -0.5; // Initial position of the first black hole
const double X2 = 0.5;  // Initial position of the second black hole

void initialize_black_holes() {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                double x = i * dx;
                double y = j * dy;
                double z = k * dz;

                double r1 = sqrt((x - X1) * (x - X1) + y * y + z * z);
                double r2 = sqrt((x - X2) * (x - X2) + y * y + z * z);

                double A_tilde_xx = 0.0, A_tilde_xy = 0.0, A_tilde_xz = 0.0;
                double A_tilde_yy = 0.0, A_tilde_yz = 0.0, A_tilde_zz = 0.0;

                if (r1 > 1e-6) {
                    double n1_x = (x - X1) / r1;
                    double n1_y = y / r1;
                    double n1_z = z / r1;
                    double P1_dot_n1 = P1 * n1_x;

                    A_tilde_xx += (3.0 / (2.0 * r1 * r1)) * (2.0 * P1_dot_n1 * n1_x - P1);
                    A_tilde_xy += (3.0 / (2.0 * r1 * r1)) * P1 * n1_y;
                    A_tilde_xz += (3.0 / (2.0 * r1 * r1)) * P1 * n1_z;
                    A_tilde_yy += (3.0 / (2.0 * r1 * r1)) * (-P1_dot_n1);
                    A_tilde_zz += (3.0 / (2.0 * r1 * r1)) * (-P1_dot_n1);
                }

                if (r2 > 1e-6) {
                    double n2_x = (x - X2) / r2;
                    double n2_y = y / r2;
                    double n2_z = z / r2;
                    double P2_dot_n2 = P2 * n2_x;

                    A_tilde_xx += (3.0 / (2.0 * r2 * r2)) * (2.0 * P2_dot_n2 * n2_x - P2);
                    A_tilde_xy += (3.0 / (2.0 * r2 * r2)) * P2 * n2_y;
                    A_tilde_xz += (3.0 / (2.0 * r2 * r2)) * P2 * n2_z;
                    A_tilde_yy += (3.0 / (2.0 * r2 * r2)) * (-P2_dot_n2);
                    A_tilde_zz += (3.0 / (2.0 * r2 * r2)) * (-P2_dot_n2);
                }
                // Assign calculated components to the grid
                A_tilde[i][j][k][0][0] = A_tilde_xx;
                A_tilde[i][j][k][0][1] = A_tilde_xy;
                A_tilde[i][j][k][0][2] = A_tilde_xz;
                A_tilde[i][j][k][1][0] = A_tilde_xy;
                A_tilde[i][j][k][1][1] = A_tilde_yy;
                A_tilde[i][j][k][1][2] = A_tilde_yz;
                A_tilde[i][j][k][2][0] = A_tilde_xz;
                A_tilde[i][j][k][2][1] = A_tilde_yz;
                A_tilde[i][j][k][2][2] = A_tilde_zz;
            }
        }
    }
}


// Tolerance for the iterative solver
const double tolerance = 1e-6;
const int max_iterations = 10000;
const double relaxation_factor = 1.0;

// Function to compute the Hamiltonian constraint and solve for phi
void solve_hamiltonian_constraint() {
    double max_residual;
    double phi_new[NX][NY][NZ];

    for (int iter = 0; iter < max_iterations; iter++) {
        max_residual = 0.0;

        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    // Compute the Laplacian of phi
                    double laplacian_phi = (phi[i+1][j][k] + phi[i-1][j][k] +
                                            phi[i][j+1][k] + phi[i][j-1][k] +
                                            phi[i][j][k+1] + phi[i][j][k-1] -
                                            6.0 * phi[i][j][k]) / (dx * dx);

                    // Compute the gradient of phi (first derivatives)
                    double dphi_dx = (phi[i+1][j][k] - phi[i-1][j][k]) / (2.0 * dx);
                    double dphi_dy = (phi[i][j+1][k] - phi[i][j-1][k]) / (2.0 * dy);
                    double dphi_dz = (phi[i][j][k+1] - phi[i][j][k-1]) / (2.0 * dz);

                    // Compute the dot product of the gradient of phi
                    double grad_phi_grad_phi = dphi_dx * dphi_dx +
                                               dphi_dy * dphi_dy +
                                               dphi_dz * dphi_dz;

                    // Placeholder for the actual Ricci scalar calculation
                    double Ricci_scalar = 0.0; // This needs to be computed from the conformal metric
                    
                    // Calculate the traceless part of the extrinsic curvature squared
                    double A_tilde_squared = 0.0;
                    for (int l = 0; l < 3; l++) {
                        for (int m = 0; m < 3; m++) {
                            A_tilde_squared += A_tilde[i][j][k][l][m] * A_tilde[i][j][k][l][m];
                        }
                    }
                    // The Hamiltonian constraint residual
                    double residual = e^{-4 * phi[i][j][k]} *
                                     (Ricci_scalar - 8.0 * laplacian_phi - 8.0 * grad_phi_grad_phi) +
                                     (2.0 / 3.0) * K[i][j][k] * K[i][j][k] - A_tilde_squared;

                    // Update the conformal factor phi
                    phi_new[i][j][k] = phi[i][j][k] - relaxation_factor * residual;

                    max_residual = std::max(max_residual, std::abs(residual));
                }
            }
        }

        // Update phi with the new values
        std::copy(&phi_new[0][0][0], &phi_new[0][0][0] + NX*NY*NZ, &phi[0][0][0]);

        // Check for convergence
        if (max_residual < tolerance) {
            std::cout << "Hamiltonian constraint converged after " << iter << " iterations." << std::endl;
            break;
        }
    }

    if (max_residual >= tolerance) {
        std::cerr << "Warning: Hamiltonian constraint did not converge!" << std::endl;
    }
}



void solve_momentum_constraint() {
    const int max_iterations = 10000;
    const double tolerance = 1e-6;

    double max_residual;

    for (int iter = 0; iter < max_iterations; iter++) {
        max_residual = 0.0;

        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    double divergence[3] = {0.0, 0.0, 0.0};

                    // Compute divergence of A_tilde
                    for (int l = 0; l < 3; l++) {
                        divergence[l] += (A_tilde[i+1][j][k][l][l] - A_tilde[i-1][j][k][l][l]) / (2.0 * dx) +
                                         (A_tilde[i][j+1][k][l][l] - A_tilde[i][j-1][k][l][l]) / (2.0 * dy) +
                                         (A_tilde[i][j][k+1][l][l] - A_tilde[i][j][k-1][l][l]) / (2.0 * dz);
                    }

                    // Compute 6 * A_tilde^ij * ∂_j φ
                    double grad_phi_correction[3] = {0.0, 0.0, 0.0};
                    double dphi_dx = (phi[i+1][j][k] - phi[i-1][j][k]) / (2.0 * dx);
                    double dphi_dy = (phi[i][j+1][k] - phi[i][j-1][k]) / (2.0 * dy);
                    double dphi_dz = (phi[i][j][k+1] - phi[i][j][k-1]) / (2.0 * dz);

                    for (int l = 0; l < 3; l++) {
                        grad_phi_correction[l] += 6.0 * (A_tilde[i][j][k][l][0] * dphi_dx +
                                                         A_tilde[i][j][k][l][1] * dphi_dy +
                                                         A_tilde[i][j][k][l][2] * dphi_dz);
                    }

                    // Compute -2/3 * ∂^i K (assuming K is constant here, so this term is zero)
                    double dK_dx = 0.0;
                    double dK_dy = 0.0;
                    double dK_dz = 0.0;
                    // Total momentum constraint residual
                    double residual[3] = {
                        divergence[0] + grad_phi_correction[0] - (2.0 / 3.0) * dK_dx,
                        divergence[1] + grad_phi_correction[1] - (2.0 / 3.0) * dK_dy,
                        divergence[2] + grad_phi_correction[2] - (2.0 / 3.0) * dK_dz
                    };

                    double total_residual = sqrt(residual[0] * residual[0] +
                                                 residual[1] * residual[1] +
                                                 residual[2] * residual[2]);

                    // Apply correction to A_tilde to reduce the residual
                    // This could be done more accurately by using a specific correction scheme
                    // Here, for simplicity, we just note the residual and could apply an adjustment

                    max_residual = std::max(max_residual, total_residual);
                }
            }
        }

        // Check for convergence
        if (max_residual < tolerance) {
            std::cout << "Momentum constraint converged after " << iter << " iterations." << std::endl;
            break;
        }
    }

    if (max_residual >= tolerance) {
        std::cerr << "Warning: Momentum constraint did not converge!" << std::endl;
    }
}

void initialize_black_holes() {
    // Initial setup of grid, black holes, etc.
    // ...
    // Solve the Hamiltonian constraint to adjust phi
    solve_hamiltonian_constraint();

    // Now solve the Momentum constraint to adjust A_tilde
    solve_momentum_constraint();
}

int ix = (i + NX) % NX;
int jx = (j + NY) % NY;
int kx = (k + NZ) % NZ;

void compute_phi_derivative(
    double phi_k[NX][NY][NZ], 
    double phi[NX][NY][NZ], 
    double alpha[NX][NY][NZ], 
    double beta_x[NX][NY][NZ], 
    double beta_y[NX][NY][NZ], 
    double beta_z[NX][NY][NZ], 
    double gamma_tilde[NX][NY][NZ][3][3], 
    double K[NX][NY][NZ], 
    double A_tilde[NX][NY][NZ][3][3], 
    double Lambda_tilde_x[NX][NY][NZ], 
    double Lambda_tilde_y[NX][NY][NZ], 
    double Lambda_tilde_z[NX][NY][NZ], 
    double B_x[NX][NY][NZ], 
    double B_y[NX][NY][NZ], 
    double B_z[NX][NY][NZ],
    double dt_factor
) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                int ip1 = (i + 1 + NX) % NX;
                int im1 = (i - 1 + NX) % NX;
                int jp1 = (j + 1 + NY) % NY;
                int jm1 = (j - 1 + NY) % NY;
                int kp1 = (k + 1 + NZ) % NZ;
                int km1 = (k - 1 + NZ) % NZ;

                double dphi_dx = (phi[ip1][j][k] - phi[im1][j][k]) / (2.0 * dx);
                double dphi_dy = (phi[i][jp1][k] - phi[i][jm1][k]) / (2.0 * dy);
                double dphi_dz = (phi[i][j][kp1] - phi[i][j][km1]) / (2.0 * dz);

                double div_beta = (beta_x[ip1][j][k] - beta_x[im1][j][k]) / (2.0 * dx) +
                                  (beta_y[i][jp1][k] - beta_y[i][jm1][k]) / (2.0 * dy) +
                                  (beta_z[i][j][kp1] - beta_z[i][j][km1]) / (2.0 * dz);

                phi_k[i][j][k] = dt_factor * (
                    (1.0 / 6.0) * div_beta 
                    - (1.0 / 6.0) * alpha[i][j][k] * K[i][j][k]
                    + beta_x[i][j][k] * dphi_dx
                    + beta_y[i][j][k] * dphi_dy
                    + beta_z[i][j][k] * dphi_dz
                );
            }
        }
    }
}



void compute_gamma_tilde_derivative(
    double gamma_tilde_k[NX][NY][NZ][3][3],
    double gamma_tilde[NX][NY][NZ][3][3],
    double alpha[NX][NY][NZ],
    double beta_x[NX][NY][NZ],
    double beta_y[NX][NY][NZ],
    double beta_z[NX][NY][NZ],
    double K[NX][NY][NZ],
    double A_tilde[NX][NY][NZ][3][3],
    double Lambda_tilde_x[NX][NY][NZ],
    double Lambda_tilde_y[NX][NY][NZ],
    double Lambda_tilde_z[NX][NY][NZ],
    double B_x[NX][NY][NZ],
    double B_y[NX][NY][NZ],
    double B_z[NX][NY][NZ],
    double dt_factor
) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        gamma_tilde_k[i][j][k][l][m] = dt_factor * (
                            -2.0 * alpha[i][j][k] * A_tilde[i][j][k][l][m]
                            + beta_x[i][j][k] * (gamma_tilde[(i+1)%NX][j][k][l][m] - gamma_tilde[(i-1+NX)%NX][j][k][l][m]) / (2.0 * dx)
                            + beta_y[i][j][k] * (gamma_tilde[i][(j+1)%NY][k][l][m] - gamma_tilde[i][(j-1+NY)%NY][k][l][m]) / (2.0 * dy)
                            + beta_z[i][j][k] * (gamma_tilde[i][j][(k+1)%NZ][l][m] - gamma_tilde[i][j][(k-1+NZ)%NZ][l][m]) / (2.0 * dz)
                        );
                    }
                }
            }
        }
    }
}



void compute_K_derivative(
    double K_k[NX][NY][NZ],
    double K[NX][NY][NZ],
    double alpha[NX][NY][NZ],
    double beta_x[NX][NY][NZ],
    double beta_y[NX][NY][NZ],
    double beta_z[NX][NY][NZ],
    double gamma_tilde[NX][NY][NZ][3][3],
    double A_tilde[NX][NY][NZ][3][3],
    double Lambda_tilde_x[NX][NY][NZ],
    double Lambda_tilde_y[NX][NY][NZ],
    double Lambda_tilde_z[NX][NY][NZ],
    double B_x[NX][NY][NZ],
    double B_y[NX][NY][NZ],
    double B_z[NX][NY][NZ],
    double dt_factor
) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                K_k[i][j][k] = dt_factor * (
                    alpha[i][j][k] * (A_tilde[i][j][k][0][0] * A_tilde[i][j][k][0][0] +
                                      A_tilde[i][j][k][1][1] * A_tilde[i][j][k][1][1] +
                                      A_tilde[i][j][k][2][2] * A_tilde[i][j][k][2][2] +
                                      2.0 * A_tilde[i][j][k][0][1] * A_tilde[i][j][k][0][1] +
                                      2.0 * A_tilde[i][j][k][0][2] * A_tilde[i][j][k][0][2] +
                                      2.0 * A_tilde[i][j][k][1][2] * A_tilde[i][j][k][1][2]) +
                    (1.0 / 3.0) * K[i][j][k] * K[i][j][k]
                );
            }
        }
    }
}



void compute_A_tilde_derivative(
    double A_tilde_k[NX][NY][NZ][3][3],
    double A_tilde[NX][NY][NZ][3][3],
    double alpha[NX][NY][NZ],
    double beta_x[NX][NY][NZ],
    double beta_y[NX][NY][NZ],
    double beta_z[NX][NY][NZ],
    double gamma_tilde[NX][NY][NZ][3][3],
    double K[NX][NY][NZ],
    double Lambda_tilde_x[NX][NY][NZ],
    double Lambda_tilde_y[NX][NY][NZ],
    double Lambda_tilde_z[NX][NY][NZ],
    double B_x[NX][NY][NZ],
    double B_y[NX][NY][NZ],
    double B_z[NX][NY][NZ],
    double dt_factor
) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        A_tilde_k[i][j][k][l][m] = dt_factor * (
                            - (2.0 / 3.0) * A_tilde[i][j][k][l][m] * K[i][j][k]
                            + alpha[i][j][k] * (K[i][j][k] * A_tilde[i][j][k][l][m]
                                                - 2.0 * A_tilde[i][j][k][l][m] * A_tilde[i][j][k][l][m])
                            + beta_x[i][j][k] * (A_tilde[(i+1)%NX][j][k][l][m] - A_tilde[(i-1+NX)%NX][j][k][l][m]) / (2.0 * dx)
                            + beta_y[i][j][k] * (A_tilde[i][(j+1)%NY][k][l][m] - A_tilde[i][(j-1+NY)%NY][k][l][m]) / (2.0 * dy)
                            + beta_z[i][j][k] * (A_tilde[i][j][(k+1)%NZ][l][m] - A_tilde[i][j][(k-1+NZ)%NZ][l][m]) / (2.0 * dz)
                        );
                    }
                }
            }
        }
    }
}



void compute_Lambda_tilde_derivative(
    double Lambda_tilde_x_k[NX][NY][NZ],
    double Lambda_tilde_y_k[NX][NY][NZ],
    double Lambda_tilde_z_k[NX][NY][NZ],
    double Lambda_tilde_x[NX][NY][NZ],
    double Lambda_tilde_y[NX][NY][NZ],
    double Lambda_tilde_z[NX][NY][NZ],
    double gamma_tilde[NX][NY][NZ][3][3],
    double alpha[NX][NY][NZ],
    double beta_x[NX][NY][NZ],
    double beta_y[NX][NY][NZ],
    double beta_z[NX][NY][NZ],
    double A_tilde[NX][NY][NZ][3][3],
    double B_x[NX][NY][NZ],
    double B_y[NX][NY][NZ],
    double B_z[NX][NY][NZ],
    double dt_factor
) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                Lambda_tilde_x_k[i][j][k] = dt_factor * (
                    alpha[i][j][k] * (
                        (A_tilde[i][j][k][0][1] * gamma_tilde[i][j][k][0][1])
                        + (A_tilde[i][j][k][0][2] * gamma_tilde[i][j][k][0][2])
                    )
                    + beta_x[i][j][k] * (Lambda_tilde_x[(i+1)%NX][j][k] - Lambda_tilde_x[(i-1+NX)%NX][j][k]) / (2.0 * dx)
                    + beta_y[i][j][k] * (Lambda_tilde_y[i][(j+1)%NY][k] - Lambda_tilde_y[i][(j-1+NY)%NY][k]) / (2.0 * dy)
                    + beta_z[i][j][k] * (Lambda_tilde_z[i][j][(k+1)%NZ] - Lambda_tilde_z[i][j][(k-1+NZ)%NZ]) / (2.0 * dz)
                );
                Lambda_tilde_y_k[i][j][k] = dt_factor * (
                    alpha[i][j][k] * (
                        (A_tilde[i][j][k][1][0] * gamma_tilde[i][j][k][1][0])
                        + (A_tilde[i][j][k][1][2] * gamma_tilde[i][j][k][1][2])
                    )
                    + beta_x[i][j][k] * (Lambda_tilde_x[(i+1)%NX][j][k] - Lambda_tilde_x[(i-1+NX)%NX][j][k]) / (2.0 * dx)
                    + beta_y[i][j][k] * (Lambda_tilde_y[i][(j+1)%NY][k] - Lambda_tilde_y[i][(j-1+NY)%NY][k]) / (2.0 * dy)
                    + beta_z[i][j][k] * (Lambda_tilde_z[i][j][(k+1)%NZ] - Lambda_tilde_z[i][j][(k-1+NZ)%NZ]) / (2.0 * dz)
                );
                Lambda_tilde_z_k[i][j][k] = dt_factor * (
                    alpha[i][j][k] * (
                        (A_tilde[i][j][k][2][0] * gamma_tilde[i][j][k][2][0])
                        + (A_tilde[i][j][k][2][1] * gamma_tilde[i][j][k][2][1])
                    )
                    + beta_x[i][j][k] * (Lambda_tilde_x[(i+1)%NX][j][k] - Lambda_tilde_x[(i-1+NX)%NX][j][k]) / (2.0 * dx)
                    + beta_y[i][j][k] * (Lambda_tilde_y[i][(j+1)%NY][k] - Lambda_tilde_y[i][(j-1+NY)%NY][k]) / (2.0 * dy)
                    + beta_z[i][j][k] * (Lambda_tilde_z[i][j][(k+1)%NZ] - Lambda_tilde_z[i][j][(k-1+NZ)%NZ]) / (2.0 * dz)
                );
            }
        }
    }
}



void compute_alpha_derivative(
    double alpha_k[NX][NY][NZ], 
    double phi[NX][NY][NZ], 
    double alpha[NX][NY][NZ], 
    double beta_x[NX][NY][NZ], 
    double beta_y[NX][NY][NZ], 
    double beta_z[NX][NY][NZ], 
    double gamma_tilde[NX][NY][NZ][3][3], 
    double K[NX][NY][NZ], 
    double A_tilde[NX][NY][NZ][3][3], 
    double Lambda_tilde_x[NX][NY][NZ], 
    double Lambda_tilde_y[NX][NY][NZ], 
    double Lambda_tilde_z[NX][NY][NZ], 
    double B_x[NX][NY][NZ], 
    double B_y[NX][NY][NZ], 
    double B_z[NX][NY][NZ],
    double dt_factor
) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                int ip1 = (i + 1 + NX) % NX;
                int im1 = (i - 1 + NX) % NX;
                int jp1 = (j + 1 + NY) % NY;
                int jm1 = (j - 1 + NY) % NY;
                int kp1 = (k + 1 + NZ) % NZ;
                int km1 = (k - 1 + NZ) % NZ;

                double dalpha_dx = (alpha[ip1][j][k] - alpha[im1][j][k]) / (2.0 * dx);
                double dalpha_dy = (alpha[i][jp1][k] - alpha[i][jm1][k]) / (2.0 * dy);
                double dalpha_dz = (alpha[i][j][kp1] - alpha[i][j][km1]) / (2.0 * dz);

                alpha_k[i][j][k] = dt_factor * (
                    -alpha[i][j][k] * K[i][j][k]
                    + beta_x[i][j][k] * dalpha_dx
                    + beta_y[i][j][k] * dalpha_dy
                    + beta_z[i][j][k] * dalpha_dz
                );
            }
        }
    }
}



void compute_beta_derivative(
    double beta_x_k[NX][NY][NZ], 
    double beta_y_k[NX][NY][NZ], 
    double beta_z_k[NX][NY][NZ], 
    double beta_x[NX][NY][NZ], 
    double beta_y[NX][NY][NZ], 
    double beta_z[NX][NY][NZ], 
    double gamma_tilde[NX][NY][NZ][3][3], 
    double K[NX][NY][NZ], 
    double A_tilde[NX][NY][NZ][3][3], 
    double Lambda_tilde_x[NX][NY][NZ], 
    double Lambda_tilde_y[NX][NY][NZ], 
    double Lambda_tilde_z[NX][NY][NZ], 
    double B_x[NX][NY][NZ], 
    double B_y[NX][NY][NZ], 
    double B_z[NX][NY][NZ],
    double dt_factor
) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                beta_x_k[i][j][k] = dt_factor * B_x[i][j][k];
                beta_y_k[i][j][k] = dt_factor * B_y[i][j][k];
                beta_z_k[i][j][k] = dt_factor * B_z[i][j][k];
            }
        }
    }
}



void compute_B_derivative(
    double B_x_k[NX][NY][NZ],
    double B_y_k[NX][NY][NZ],
    double B_z_k[NX][NY][NZ],
    double B_x[NX][NY][NZ],
    double B_y[NX][NY][NZ],
    double B_z[NX][NY][NZ],
    double gamma_tilde[NX][NY][NZ][3][3],
    double alpha[NX][NY][NZ],
    double beta_x[NX][NY][NZ],
    double beta_y[NX][NY][NZ],
    double beta_z[NX][NY][NZ],
    double Lambda_tilde_x[NX][NY][NZ],
    double Lambda_tilde_y[NX][NY][NZ],
    double Lambda_tilde_z[NX][NY][NZ],
    double dt_factor
) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                B_x_k[i][j][k] = dt_factor * (
                    alpha[i][j][k] * (Lambda_tilde_x[i][j][k] - beta_x[i][j][k])
                );
                B_y_k[i][j][k] = dt_factor * (
                    alpha[i][j][k] * (Lambda_tilde_y[i][j][k] - beta_y[i][j][k])
                );
                B_z_k[i][j][k] = dt_factor * (
                    alpha[i][j][k] * (Lambda_tilde_z[i][j][k] - beta_z[i][j][k])
                );
            }
        }
    }
}







void rk4_step(double dt) {
    // Temporary arrays for k1, k2, k3, k4 steps for all variables
    double phi_k1[NX][NY][NZ], phi_k2[NX][NY][NZ], phi_k3[NX][NY][NZ], phi_k4[NX][NY][NZ];
    double alpha_k1[NX][NY][NZ], alpha_k2[NX][NY][NZ], alpha_k3[NX][NY][NZ], alpha_k4[NX][NY][NZ];
    double beta_x_k1[NX][NY][NZ], beta_x_k2[NX][NY][NZ], beta_x_k3[NX][NY][NZ], beta_x_k4[NX][NY][NZ];
    double beta_y_k1[NX][NY][NZ], beta_y_k2[NX][NY][NZ], beta_y_k3[NX][NY][NZ], beta_y_k4[NX][NY][NZ];
    double beta_z_k1[NX][NY][NZ], beta_z_k2[NX][NY][NZ], beta_z_k3[NX][NY][NZ], beta_z_k4[NX][NY][NZ];
    double gamma_tilde_k1[NX][NY][NZ][3][3], gamma_tilde_k2[NX][NY][NZ][3][3], gamma_tilde_k3[NX][NY][NZ][3][3], gamma_tilde_k4[NX][NY][NZ][3][3];
    double K_k1[NX][NY][NZ], K_k2[NX][NY][NZ], K_k3[NX][NY][NZ], K_k4[NX][NY][NZ];
    double A_tilde_k1[NX][NY][NZ][3][3], A_tilde_k2[NX][NY][NZ][3][3], A_tilde_k3[NX][NY][NZ][3][3], A_tilde_k4[NX][NY][NZ][3][3];
    double Lambda_tilde_x_k1[NX][NY][NZ], Lambda_tilde_x_k2[NX][NY][NZ], Lambda_tilde_x_k3[NX][NY][NZ], Lambda_tilde_x_k4[NX][NY][NZ];
    double Lambda_tilde_y_k1[NX][NY][NZ], Lambda_tilde_y_k2[NX][NY][NZ], Lambda_tilde_y_k3[NX][NY][NZ], Lambda_tilde_y_k4[NX][NY][NZ];
    double Lambda_tilde_z_k1[NX][NY][NZ], Lambda_tilde_z_k2[NX][NY][NZ], Lambda_tilde_z_k3[NX][NY][NZ], Lambda_tilde_z_k4[NX][NY][NZ];
    double B_x_k1[NX][NY][NZ], B_x_k2[NX][NY][NZ], B_x_k3[NX][NY][NZ], B_x_k4[NX][NY][NZ];
    double B_y_k1[NX][NY][NZ], B_y_k2[NX][NY][NZ], B_y_k3[NX][NY][NZ], B_y_k4[NX][NY][NZ];
    double B_z_k1[NX][NY][NZ], B_z_k2[NX][NY][NZ], B_z_k3[NX][NY][NZ], B_z_k4[NX][NY][NZ];

    // Temporary arrays to store intermediate values
    double phi_temp[NX][NY][NZ], alpha_temp[NX][NY][NZ], beta_x_temp[NX][NY][NZ], beta_y_temp[NX][NY][NZ], beta_z_temp[NX][NY][NZ];
    double gamma_tilde_temp[NX][NY][NZ][3][3], K_temp[NX][NY][NZ], A_tilde_temp[NX][NY][NZ][3][3];
    double Lambda_tilde_x_temp[NX][NY][NZ], Lambda_tilde_y_temp[NX][NY][NZ], Lambda_tilde_z_temp[NX][NY][NZ];
    double B_x_temp[NX][NY][NZ], B_y_temp[NX][NY][NZ], B_z_temp[NX][NY][NZ];

    // Step 1: Compute k1 for all variables using the original variables
    compute_phi_derivative(phi_k1, phi, alpha, beta_x, beta_y, beta_z, gamma_tilde, K, A_tilde, Lambda_tilde_x, Lambda_tilde_y, Lambda_tilde_z, B_x, B_y, B_z, dt);
    compute_alpha_derivative(alpha_k1, phi, alpha, beta_x, beta_y, beta_z, gamma_tilde, K, A_tilde, Lambda_tilde_x, Lambda_tilde_y, Lambda_tilde_z, B_x, B_y, B_z, dt);
    compute_beta_derivative(beta_x_k1, beta_y_k1, beta_z_k1, beta_x, beta_y, beta_z, gamma_tilde, K, A_tilde, Lambda_tilde_x, Lambda_tilde_y, Lambda_tilde_z, B_x, B_y, B_z, dt);
    compute_gamma_tilde_derivative(gamma_tilde_k1, gamma_tilde, alpha, beta_x, beta_y, beta_z, K, A_tilde, Lambda_tilde_x, Lambda_tilde_y, Lambda_tilde_z, B_x, B_y, B_z, dt);
    compute_K_derivative(K_k1, K, alpha, beta_x, beta_y, beta_z, gamma_tilde, A_tilde, Lambda_tilde_x, Lambda_tilde_y, Lambda_tilde_z, B_x, B_y, B_z, dt);
    compute_A_tilde_derivative(A_tilde_k1, A_tilde, alpha, beta_x, beta_y, beta_z, gamma_tilde, K, Lambda_tilde_x, Lambda_tilde_y, Lambda_tilde_z, B_x, B_y, B_z, dt);
    compute_Lambda_tilde_derivative(Lambda_tilde_x_k1, Lambda_tilde_y_k1, Lambda_tilde_z_k1, Lambda_tilde_x, Lambda_tilde_y, Lambda_tilde_z, gamma_tilde, alpha, beta_x, beta_y, beta_z, A_tilde, B_x, B_y, B_z, dt);
    compute_B_derivative(B_x_k1, B_y_k1, B_z_k1, B_x, B_y, B_z, gamma_tilde, alpha, beta_x, beta_y, beta_z, Lambda_tilde_x, Lambda_tilde_y, Lambda_tilde_z, dt);

    // Step 2: Compute k2 using k1 and the original variables
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                phi_temp[i][j][k] = phi[i][j][k] + 0.5 * phi_k1[i][j][k];
                alpha_temp[i][j][k] = alpha[i][j][k] + 0.5 * alpha_k1[i][j][k];
                beta_x_temp[i][j][k] = beta_x[i][j][k] + 0.5 * beta_x_k1[i][j][k];
                beta_y_temp[i][j][k] = beta_y[i][j][k] + 0.5 * beta_y_k1[i][j][k];
                beta_z_temp[i][j][k] = beta_z[i][j][k] + 0.5 * beta_z_k1[i][j][k];
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        gamma_tilde_temp[i][j][k][l][m] = gamma_tilde[i][j][k][l][m] + 0.5 * dt * gamma_tilde_k1[i][j][k][l][m];
                        A_tilde_temp[i][j][k][l][m] = A_tilde[i][j][k][l][m] + 0.5 * dt * A_tilde_k1[i][j][k][l][m];
                    }
                }
                K_temp[i][j][k] = K[i][j][k] + 0.5 * dt * K_k1[i][j][k];
                Lambda_tilde_x_temp[i][j][k] = Lambda_tilde_x[i][j][k] + 0.5 * dt * Lambda_tilde_x_k1[i][j][k];
                Lambda_tilde_y_temp[i][j][k] = Lambda_tilde_y[i][j][k] + 0.5 * dt * Lambda_tilde_y_k1[i][j][k];
                Lambda_tilde_z_temp[i][j][k] = Lambda_tilde_z[i][j][k] + 0.5 * dt * Lambda_tilde_z_k1[i][j][k];
                B_x_temp[i][j][k] = B_x[i][j][k] + 0.5 * dt * B_x_k1[i][j][k];
                B_y_temp[i][j][k] = B_y[i][j][k] + 0.5 * dt * B_y_k1[i][j][k];
                B_z_temp[i][j][k] = B_z[i][j][k] + 0.5 * dt * B_z_k1[i][j][k];
            }
        }
    }
    compute_phi_derivative(phi_k2, phi_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, K_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_alpha_derivative(alpha_k2, phi_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, K_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_beta_derivative(beta_x_k2, beta_y_k2, beta_z_k2, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, K_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_gamma_tilde_derivative(gamma_tilde_k2, gamma_tilde_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, K_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_K_derivative(K_k2, K_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_A_tilde_derivative(A_tilde_k2, A_tilde_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, K_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_Lambda_tilde_derivative(Lambda_tilde_x_k2, Lambda_tilde_y_k2, Lambda_tilde_z_k2, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, gamma_tilde_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, A_tilde_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_B_derivative(B_x_k2, B_y_k2, B_z_k2, B_x_temp, B_y_temp, B_z_temp, gamma_tilde_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, dt);

    // Step 3: Compute k3 using k2 and the original variables
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                phi_temp[i][j][k] = phi[i][j][k] + 0.5 * phi_k2[i][j][k];
                alpha_temp[i][j][k] = alpha[i][j][k] + 0.5 * alpha_k2[i][j][k];
                beta_x_temp[i][j][k] = beta_x[i][j][k] + 0.5 * beta_x_k2[i][j][k];
                beta_y_temp[i][j][k] = beta_y[i][j][k] + 0.5 * beta_y_k2[i][j][k];
                beta_z_temp[i][j][k] = beta_z[i][j][k] + 0.5 * beta_z_k2[i][j][k];
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        gamma_tilde_temp[i][j][k][l][m] = gamma_tilde[i][j][k][l][m] + 0.5 * dt * gamma_tilde_k2[i][j][k][l][m];
                        A_tilde_temp[i][j][k][l][m] = A_tilde[i][j][k][l][m] + 0.5 * dt * A_tilde_k2[i][j][k][l][m];
                    }
                }
                K_temp[i][j][k] = K[i][j][k] + 0.5 * dt * K_k2[i][j][k];
                Lambda_tilde_x_temp[i][j][k] = Lambda_tilde_x[i][j][k] + 0.5 * dt * Lambda_tilde_x_k2[i][j][k];
                Lambda_tilde_y_temp[i][j][k] = Lambda_tilde_y[i][j][k] + 0.5 * dt * Lambda_tilde_y_k2[i][j][k];
                Lambda_tilde_z_temp[i][j][k] = Lambda_tilde_z[i][j][k] + 0.5 * dt * Lambda_tilde_z_k2[i][j][k];
                B_x_temp[i][j][k] = B_x[i][j][k] + 0.5 * dt * B_x_k2[i][j][k];
                B_y_temp[i][j][k] = B_y[i][j][k] + 0.5 * dt * B_y_k2[i][j][k];
                B_z_temp[i][j][k] = B_z[i][j][k] + 0.5 * dt * B_z_k2[i][j][k];
            }
        }
    }
    compute_phi_derivative(phi_k3, phi_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, K_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_alpha_derivative(alpha_k3, phi_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, K_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_beta_derivative(beta_x_k3, beta_y_k3, beta_z_k3, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, K_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_gamma_tilde_derivative(gamma_tilde_k3, gamma_tilde_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, K_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_K_derivative(K_k3, K_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_A_tilde_derivative(A_tilde_k3, A_tilde_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, K_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_Lambda_tilde_derivative(Lambda_tilde_x_k3, Lambda_tilde_y_k3, Lambda_tilde_z_k3, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, gamma_tilde_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, A_tilde_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_B_derivative(B_x_k3, B_y_k3, B_z_k3, B_x_temp, B_y_temp, B_z_temp, gamma_tilde_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, dt);

    // Step 4: Compute k4 using k3 and the original variables
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                phi_temp[i][j][k] = phi[i][j][k] + phi_k3[i][j][k];
                alpha_temp[i][j][k] = alpha[i][j][k] + alpha_k3[i][j][k];
                beta_x_temp[i][j][k] = beta_x[i][j][k] + beta_x_k3[i][j][k];
                beta_y_temp[i][j][k] = beta_y[i][j][k] + beta_y_k3[i][j][k];
                beta_z_temp[i][j][k] = beta_z[i][j][k] + beta_z_k3[i][j][k];
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        gamma_tilde_temp[i][j][k][l][m] = gamma_tilde[i][j][k][l][m] + dt * gamma_tilde_k3[i][j][k][l][m];
                        A_tilde_temp[i][j][k][l][m] = A_tilde[i][j][k][l][m] + dt * A_tilde_k3[i][j][k][l][m];
                    }
                }
                K_temp[i][j][k] = K[i][j][k] + dt * K_k3[i][j][k];
                Lambda_tilde_x_temp[i][j][k] = Lambda_tilde_x[i][j][k] + dt * Lambda_tilde_x_k3[i][j][k];
                Lambda_tilde_y_temp[i][j][k] = Lambda_tilde_y[i][j][k] + dt * Lambda_tilde_y_k3[i][j][k];
                Lambda_tilde_z_temp[i][j][k] = Lambda_tilde_z[i][j][k] + dt * Lambda_tilde_z_k3[i][j][k];
                B_x_temp[i][j][k] = B_x[i][j][k] + dt * B_x_k3[i][j][k];
                B_y_temp[i][j][k] = B_y[i][j][k] + dt * B_y_k3[i][j][k];
                B_z_temp[i][j][k] = B_z[i][j][k] + dt * B_z_k3[i][j][k];
            }
        }
    }
    compute_phi_derivative(phi_k4, phi_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, K_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_alpha_derivative(alpha_k4, phi_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, K_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_beta_derivative(beta_x_k4, beta_y_k4, beta_z_k4, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, K_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_gamma_tilde_derivative(gamma_tilde_k4, gamma_tilde_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, K_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_K_derivative(K_k4, K_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, A_tilde_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_A_tilde_derivative(A_tilde_k4, A_tilde_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, gamma_tilde_temp, K_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_Lambda_tilde_derivative(Lambda_tilde_x_k4, Lambda_tilde_y_k4, Lambda_tilde_z_k4, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, gamma_tilde_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, A_tilde_temp, B_x_temp, B_y_temp, B_z_temp, dt);
    compute_B_derivative(B_x_k4, B_y_k4, B_z_k4, B_x_temp, B_y_temp, B_z_temp, gamma_tilde_temp, alpha_temp, beta_x_temp, beta_y_temp, beta_z_temp, Lambda_tilde_x_temp, Lambda_tilde_y_temp, Lambda_tilde_z_temp, dt);

    // Final update step using all k1, k2, k3, k4
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                phi[i][j][k] += dt / 6.0 * (phi_k1[i][j][k] + 2.0 * phi_k2[i][j][k] + 2.0 * phi_k3[i][j][k] + phi_k4[i][j][k]);
                alpha[i][j][k] += dt / 6.0 * (alpha_k1[i][j][k] + 2.0 * alpha_k2[i][j][k] + 2.0 * alpha_k3[i][j][k] + alpha_k4[i][j][k]);
                beta_x[i][j][k] += dt / 6.0 * (beta_x_k1[i][j][k] + 2.0 * beta_x_k2[i][j][k] + 2.0 * beta_x_k3[i][j][k] + beta_x_k4[i][j][k]);
                beta_y[i][j][k] += dt / 6.0 * (beta_y_k1[i][j][k] + 2.0 * beta_y_k2[i][j][k] + 2.0 * beta_y_k3[i][j][k] + beta_y_k4[i][j][k]);
                beta_z[i][j][k] += dt / 6.0 * (beta_z_k1[i][j][k] + 2.0 * beta_z_k2[i][j][k] + 2.0 * beta_z_k3[i][j][k] + beta_z_k4[i][j][k]);
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        gamma_tilde[i][j][k][l][m] += dt / 6.0 * (gamma_tilde_k1[i][j][k][l][m] + 2.0 * gamma_tilde_k2[i][j][k][l][m] + 2.0 * gamma_tilde_k3[i][j][k][l][m] + gamma_tilde_k4[i][j][k][l][m]);
                        A_tilde[i][j][k][l][m] += dt / 6.0 * (A_tilde_k1[i][j][k][l][m] + 2.0 * A_tilde_k2[i][j][k][l][m] + 2.0 * A_tilde_k3[i][j][k][l][m] + A_tilde_k4[i][j][k][l][m]);
                    }
                }
                K[i][j][k] += dt / 6.0 * (K_k1[i][j][k] + 2.0 * K_k2[i][j][k] + 2.0 * K_k3[i][j][k] + K_k4[i][j][k]);
                Lambda_tilde_x[i][j][k] += dt / 6.0 * (Lambda_tilde_x_k1[i][j][k] + 2.0 * Lambda_tilde_x_k2[i][j][k] + 2.0 * Lambda_tilde_x_k3[i][j][k] + Lambda_tilde_x_k4[i][j][k]);
                Lambda_tilde_y[i][j][k] += dt / 6.0 * (Lambda_tilde_y_k1[i][j][k] + 2.0 * Lambda_tilde_y_k2[i][j][k] + 2.0 * Lambda_tilde_y_k3[i][j][k] + Lambda_tilde_y_k4[i][j][k]);
                Lambda_tilde_z[i][j][k] += dt / 6.0 * (Lambda_tilde_z_k1[i][j][k] + 2.0 * Lambda_tilde_z_k2[i][j][k] + 2.0 * Lambda_tilde_z_k3[i][j][k] + Lambda_tilde_z_k4[i][j][k]);
                B_x[i][j][k] += dt / 6.0 * (B_x_k1[i][j][k] + 2.0 * B_x_k2[i][j][k] + 2.0 * B_x_k3[i][j][k] + B_x_k4[i][j][k]);
                B_y[i][j][k] += dt / 6.0 * (B_y_k1[i][j][k] + 2.0 * B_y_k2[i][j][k] + 2.0 * B_y_k3[i][j][k] + B_y_k4[i][j][k]);
                B_z[i][j][k] += dt / 6.0 * (B_z_k1[i][j][k] + 2.0 * B_z_k2[i][j][k] + 2.0 * B_z_k3[i][j][k] + B_z_k4[i][j][k]);
            }
        }
    }
}






