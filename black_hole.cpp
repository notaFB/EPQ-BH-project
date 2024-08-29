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

void evolve_phi(double dt) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                // Wrap-around indices
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

                double dphi_dt = (1.0 / 6.0) * div_beta - (1.0 / 6.0) * alpha[i][j][k] * K[i][j][k]
                                 + beta_x[i][j][k] * dphi_dx
                                 + beta_y[i][j][k] * dphi_dy
                                 + beta_z[i][j][k] * dphi_dz;

                phi[i][j][k] += dt * dphi_dt;
            }
        }
    }
}

void evolve_gamma_tilde(double dt) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                int ip1 = (i + 1 + NX) % NX;
                int im1 = (i - 1 + NX) % NX;
                int jp1 = (j + 1 + NY) % NY;
                int jm1 = (j - 1 + NY) % NY;
                int kp1 = (k + 1 + NZ) % NZ;
                int km1 = (k - 1 + NZ) % NZ;

                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        double dbeta_dx = (beta_x[ip1][j][k] - beta_x[im1][j][k]) / (2.0 * dx);
                        double dbeta_dy = (beta_y[i][jp1][k] - beta_y[i][jm1][k]) / (2.0 * dy);
                        double dbeta_dz = (beta_z[i][j][kp1] - beta_z[i][j][km1]) / (2.0 * dz);

                        double div_beta = dbeta_dx + dbeta_dy + dbeta_dz;

                        double dgamma_tilde_dt = -2.0 * alpha[i][j][k] * A_tilde[i][j][k][l][m]
                            - (2.0 / 3.0) * gamma_tilde[i][j][k][l][m] * div_beta
                            + beta_x[i][j][k] * (gamma_tilde[ip1][j][k][l][m] - gamma_tilde[im1][j][k][l][m]) / (2.0 * dx)
                            + beta_y[i][j][k] * (gamma_tilde[i][jp1][k][l][m] - gamma_tilde[i][jm1][k][l][m]) / (2.0 * dy)
                            + beta_z[i][j][k] * (gamma_tilde[i][j][kp1][l][m] - gamma_tilde[i][j][km1][l][m]) / (2.0 * dz);

                        gamma_tilde[i][j][k][l][m] += dt * dgamma_tilde_dt;
                    }
                }
            }
        }
    }
}

void evolve_K(double dt) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                int ip1 = (i + 1 + NX) % NX;
                int im1 = (i - 1 + NX) % NX;
                int jp1 = (j + 1 + NY) % NY;
                int jm1 = (j - 1 + NY) % NY;
                int kp1 = (k + 1 + NZ) % NZ;
                int km1 = (k - 1 + NZ) % NZ;

                double laplacian_alpha = (alpha[ip1][j][k] - 2.0 * alpha[i][j][k] + alpha[im1][j][k]) / (dx * dx) +
                                         (alpha[i][jp1][k] - 2.0 * alpha[i][j][k] + alpha[i][jm1][k]) / (dy * dy) +
                                         (alpha[i][j][kp1] - 2.0 * alpha[i][j][k] + alpha[i][j][km1]) / (dz * dz);

                double A_tilde_squared = 0.0;
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        A_tilde_squared += A_tilde[i][j][k][l][m] * A_tilde[i][j][k][l][m];
                    }
                }

                double dK_dt = alpha[i][j][k] * (A_tilde_squared + (1.0 / 3.0) * K[i][j][k] * K[i][j][k])
                    - laplacian_alpha
                    + beta_x[i][j][k] * (K[ip1][j][k] - K[im1][j][k]) / (2.0 * dx)
                    + beta_y[i][j][k] * (K[i][jp1][k] - K[i][jm1][k]) / (2.0 * dy)
                    + beta_z[i][j][k] * (K[i][j][kp1] - K[i][j][km1]) / (2.0 * dz);

                K[i][j][k] += dt * dK_dt;
            }
        }
    }
}

void evolve_A_tilde(double dt) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                int ip1 = (i + 1 + NX) % NX;
                int im1 = (i - 1 + NX) % NX;
                int jp1 = (j + 1 + NY) % NY;
                int jm1 = (j - 1 + NY) % NY;
                int kp1 = (k + 1 + NZ) % NZ;
                int km1 = (k - 1 + NZ) % NZ;

                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        double dbeta_dx = (beta_x[ip1][j][k] - beta_x[im1][j][k]) / (2.0 * dx);
                        double dbeta_dy = (beta_y[i][jp1][k] - beta_y[i][jm1][k]) / (2.0 * dy);
                        double dbeta_dz = (beta_z[i][j][kp1] - beta_z[i][j][km1]) / (2.0 * dz);

                        double div_beta = dbeta_dx + dbeta_dy + dbeta_dz;

                        double dA_tilde_dt = - (2.0 / 3.0) * A_tilde[i][j][k][l][m] * div_beta
                            + alpha[i][j][k] * (
                                K[i][j][k] * A_tilde[i][j][k][l][m]
                                - 2.0 * (A_tilde[i][j][k][l][0] * A_tilde[i][j][k][0][m])
                            )
                            + 4.0 * exp(-4.0 * phi[i][j][k]) * (
                                alpha[i][j][k] * Ricci_ij[i][j][k][l][m]
                                - (alpha[ip1][j][k] - 2.0 * alpha[i][j][k] + alpha[im1][j][k]) / (dx * dx)
                            )
                            + beta_x[i][j][k] * (A_tilde[ip1][j][k][l][m] - A_tilde[im1][j][k][l][m]) / (2.0 * dx)
                            + beta_y[i][j][k] * (A_tilde[i][jp1][k][l][m] - A_tilde[i][jm1][k][l][m]) / (2.0 * dy)
                            + beta_z[i][j][k] * (A_tilde[i][j][kp1][l][m] - A_tilde[i][j][km1][l][m]) / (2.0 * dz);

                        A_tilde[i][j][k][l][m] += dt * dA_tilde_dt;
                    }
                }
            }
        }
    }
}

void evolve_Lambda_tilde(double dt) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                int ip1 = (i + 1 + NX) % NX;
                int im1 = (i - 1 + NX) % NX;
                int jp1 = (j + 1 + NY) % NY;
                int jm1 = (j - 1 + NY) % NY;
                int kp1 = (k + 1 + NZ) % NZ;
                int km1 = (k - 1 + NZ) % NZ;

                double dGamma_tilde_dt = (
                    (beta_x[ip1][j][k] - 2.0 * beta_x[i][j][k] + beta_x[im1][j][k]) / (dx * dx) +
                    (beta_y[i][jp1][k] - 2.0 * beta_y[i][j][k] + beta_y[i][jm1][k]) / (dy * dy) +
                    (beta_z[i][j][kp1] - 2.0 * beta_z[i][j][k] + beta_z[i][j][km1]) / (dz * dz)
                ) - 2.0 * A_tilde[i][j][k][0][0] * alpha[i][j][k]
                + beta_x[i][j][k] * (Gamma_tilde_x[ip1][j][k] - Gamma_tilde_x[im1][j][k]) / (2.0 * dx)
                + beta_y[i][j][k] * (Gamma_tilde_y[i][jp1][k] - Gamma_tilde_y[i][jm1][k]) / (2.0 * dy)
                + beta_z[i][j][k] * (Gamma_tilde_z[i][j][kp1] - Gamma_tilde_z[i][j][km1]) / (2.0 * dz);

                Lambda_tilde_x[i][j][k] += dt * dGamma_tilde_dt;
                Lambda_tilde_y[i][j][k] += dt * dGamma_tilde_dt;
                Lambda_tilde_z[i][j][k] += dt * dGamma_tilde_dt;
            }
        }
    }
}

void evolve_alpha(double dt) {
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

                double dalpha_dt = -alpha[i][j][k] * alpha[i][j][k] * K[i][j][k]
                                   + beta_x[i][j][k] * dalpha_dx
                                   + beta_y[i][j][k] * dalpha_dy
                                   + beta_z[i][j][k] * dalpha_dz;

                alpha[i][j][k] += dt * dalpha_dt;
            }
        }
    }
}

void evolve_beta(double dt) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                beta_x[i][j][k] += dt * B_x[i][j][k];
                beta_y[i][j][k] += dt * B_y[i][j][k];
                beta_z[i][j][k] += dt * B_z[i][j][k];
            }
        }
    }
}

void evolve_B(double dt, double eta) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                int ip1 = (i + 1 + NX) % NX;
                int im1 = (i - 1 + NX) % NX;
                int jp1 = (j + 1 + NY) % NY;
                int jm1 = (j - 1 + NY) % NY;
                int kp1 = (k + 1 + NZ) % NZ;
                int km1 = (k - 1 + NZ) % NZ;

                double dGamma_tilde_dx = (Gamma_tilde_x[ip1][j][k] - Gamma_tilde_x[im1][j][k]) / (2.0 * dx);
                double dGamma_tilde_dy = (Gamma_tilde_y[i][jp1][k] - Gamma_tilde_y[i][jm1][k]) / (2.0 * dy);
                double dGamma_tilde_dz = (Gamma_tilde_z[i][j][kp1] - Gamma_tilde_z[i][j][km1]) / (2.0 * dz);

                B_x[i][j][k] += dt * (dGamma_tilde_dx - eta * B_x[i][j][k] + beta_x[i][j][k] * dGamma_tilde_dx);
                B_y[i][j][k] += dt * (dGamma_tilde_dy - eta * B_y[i][j][k] + beta_y[i][j][k] * dGamma_tilde_dy);
                B_z[i][j][k] += dt * (dGamma_tilde_dz - eta * B_z[i][j][k] + beta_z[i][j][k] * dGamma_tilde_dz);
            }
        }
    }
}

void evolve_bssn(double dt, double eta) {
    evolve_phi(dt);
    evolve_gamma_tilde(dt);
    evolve_K(dt);
    evolve_A_tilde(dt);
    evolve_Lambda_tilde(dt);
    evolve_alpha(dt);
    evolve_beta(dt);
    evolve_B(dt, eta);
}


