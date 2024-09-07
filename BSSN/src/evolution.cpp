
#include"evolution.hpp"
#include "utilities.hpp"
#include "geometry.hpp"
#include <cmath>

void Evolution::RK4_step()
{

    //Apply this (getting the data from FieldData):
    /*
        # Initial condition: v1 contains the state vector s

        # Step 1: Compute k1 and store it in v2
        v2 = f(v1)           # v2 = f(s) (k1 = f(s))

        # Step 2: Compute k2 using a temporary update to v3, store the result in v3
        v3 = v1              # Copy v1 (s) into v3 to preserve original v1
        v2 *= 0.5            # Scale k1 by 0.5 (v2 = 0.5 * k1)
        v3 += v2             # v3 = s + 0.5 * k1
        v2 = f(v3)           # v2 = f(s + 0.5 * k1) (k2)

        # Step 3: Compute k3 using another temporary update to v4, store the result in v4
        v4 = v1              # Copy v1 (s) into v4 to preserve original v1
        v2 *= 0.5            # Scale k2 by 0.5 (v2 = 0.5 * k2)
        v4 += v2             # v4 = s + 0.5 * k2
        v2 = f(v4)           # v2 = f(s + 0.5 * k2) (k3)

        # Step 4: Compute k4 using v3 (reused), store the result in v3
        v3 = v1              # Copy v1 (s) into v3 again
        v2 *= 1.0            # Scale k3 by 1.0 (v2 = k3)
        v3 += v2             # v3 = s + k3
        v2 = f(v3)           # v2 = f(s + k3) (k4)

        # Step 5: Final combination of k1, k2, k3, k4 to update v1
        v3 = v1              # Copy v1 (s) into v3 for final updates
        v2 = v2              # v2 is already k4, keep it unchanged
        v3 += v2             # v3 += k4
        v2 = v2              # Copy k1 (stored earlier in v2) again for the combination
        v2 += v2             # Scale k1 by 2 and add it
        v4

    */

    //Get the data from FieldData
    FieldData& fieldData = FieldData::getInstance();
    SpatialSlice& v1 = fieldData.getV1Slice(); // v1 stores the current state vector s
    SpatialSlice& v2 = fieldData.getV2Slice(); // v2 will be used for k1, k2, k3, k4
    SpatialSlice& v3 = fieldData.getV3Slice(); // v3 will be used as temporary storage
    SpatialSlice& v4 = fieldData.getV4Slice(); // v4 will be used as temporary storage

    double halfStep = 0.5;
    double fullStep = 1.0;

    // Step 1: Compute k1 and store it in v2
    F(v1, v2);  // v2 = f(v1), which is k1
    Geometry::compute(v2);

    // Step 2: Compute k2 using v3 as temporary storage
    v3 = v1;    // Copy v1 to v3
    v2 *= halfStep; // Scale k1 by 0.5 (v2 = 0.5 * k1)
    v3 += v2;  // v3 = s + 0.5 * k1
    F(v3, v2); // v2 = f(s + 0.5 * k1), which is k2
    Geometry::compute(v2);

    // Step 3: Compute k3 using v4 as temporary storage
    v4 = v1;    // Copy v1 to v4
    v2 *= halfStep; // Scale k2 by 0.5 (v2 = 0.5 * k2)
    v4 += v2;  // v4 = s + 0.5 * k2
    F(v4, v2); // v2 = f(s + 0.5 * k2), which is k3
    Geometry::compute(v2);

    // Step 4: Compute k4 using v3 (reused)
    v3 = v1;    // Copy v1 to v3
    v2 *= fullStep; // Scale k3 by 1.0 (v2 = k3)
    v3 += v2;  // v3 = s + k3
    F(v3, v2); // v2 = f(s + k3), which is k4
    Geometry::compute(v2);

    // Step 5: Final combination of k1, k2, k3, k4 to update v1
    v3 = v2;    // v3 = k4 (copy k4 from v2)
    v3 *= fullStep; // No scaling needed for k4, keep it as is

    // Now, use v4 to compute the final sum
    v4 = v1;    // Start with v1 (the original state)
    v2 = fieldData.getV2Slice(); // Get k1 again
    v2 *= 2.0;  // Scale k1 by 2.0
    v4 += v2;  // Add 2 * k1 to v4

    v3 += v4;  // Add k4 to v3 (final result)

    v1 = v3;   // Update v1 with the final result of the RK4 step

    // REVISE THIS

}

void Evolution::F_phi(SpatialSlice& slice, SpatialSlice& sliceOut) {
    // dphi/dt = -alpha*K/6

    ScalarField& phi_out = sliceOut.phi;

    ScalarField& alpha = slice.alpha;
    ScalarField& K = slice.K;
    

    for (int nx = 0; nx < N; nx++) {
        for (int ny = 0; ny < N; ny++) {
            for (int nz = 0; nz < N; nz++) {
                phi_out(nx, ny, nz) = -alpha(nx, ny, nz) * K(nx, ny, nz) / 6.0;
            }
        }
    }


}

void Evolution::F_K(SpatialSlice& slice, SpatialSlice& sliceOut) {
    // \partial_t K = -\gamma^{ij} \left( \partial_j \partial_i \alpha - \Gamma^k_{ji} \partial_k \alpha \right) + \alpha \left( \tilde{A}_{ij} \tilde{A}^{ij} + \frac{1}{3} K^2 \right)

    ScalarField& K_out = sliceOut.K;

    ScalarField& alpha = slice.alpha;
    TensorField& gammaInv = slice.gammaInv;
    Tensor3Field& Christoffel = slice.Christoffel;
    TensorField& ATilde = slice.ATilde;
    TensorField& ATildeInv = slice.ATildeInv;
    ScalarField& K = slice.K;

    for (int nx = 0; nx < N; nx++) {
        for (int ny = 0; ny < N; ny++) {
            for (int nz = 0; nz < N; nz++) {
                
                K_out(nx, ny, nz) = 0.0;

                //First term: gamma^{ij} \left( \partial_j \partial_i \alpha - \Gamma^k_{ji} \partial_k \alpha \right)
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        for(int k = 0; k < 3; k++) {
                            K_out(nx, ny, nz) += gammaInv(nx, ny, nz, i, j) * (SECOND_PARTIAL_SCALAR(alpha, i, j) - Christoffel(nx, ny, nz, k, j, i) * PARTIAL_SCALAR(alpha, k));
                        }
                    }
                }

                //Second term: \alpha \left( \tilde{A}_{ij} \tilde{A}^{ij} + \frac{1}{3} K^2 \right)
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        K_out(nx, ny, nz) += alpha(nx, ny, nz) * (ATilde(nx, ny, nz, i, j) * ATildeInv(nx, ny, nz, i, j) + 1.0/3.0 * K(nx, ny, nz) * K(nx, ny, nz));
                    }
                }

            }
        }
    }

    
}

void Evolution::F_gammaTilde(SpatialSlice& slice, SpatialSlice& sliceOut) {
    // \partial_t \tilde{\gamma}_{ij} = -2 \alpha \tilde{A}_{ij}

    TensorField& gammaTilde_out = sliceOut.gammaTilde;

    ScalarField& alpha = slice.alpha;
    TensorField& ATilde = slice.ATilde;

    for (int nx = 0; nx < N; nx++) {
        for (int ny = 0; ny < N; ny++) {
            for (int nz = 0; nz < N; nz++) {
                
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        gammaTilde_out(nx, ny, nz, i, j) = -2.0 * alpha(nx, ny, nz) * ATilde(nx, ny, nz, i, j);
                    }
                }

            }
        }
    }
}

void Evolution::F_alpha(SpatialSlice& slice, SpatialSlice& sliceOut) {
    // \partial_t \alpha = -2 \alpha K

    ScalarField& alpha_out = sliceOut.alpha;

    ScalarField& alpha = slice.alpha;
    ScalarField& K = slice.K;

    for (int nx = 0; nx < N; nx++) {
        for (int ny = 0; ny < N; ny++) {
            for (int nz = 0; nz < N; nz++) {

                alpha_out(nx, ny, nz) = -2.0 * alpha(nx, ny, nz) * K(nx, ny, nz);
                
            }
        }
    }
}

void Evolution::F_beta(SpatialSlice& slice, SpatialSlice& sliceOut) {
    // \partial_t \beta^i = B^i

    VectorField& beta_out = sliceOut.beta;

    VectorField& B = slice.B;

    for (int nx = 0; nx < N; nx++) {
        for (int ny = 0; ny < N; ny++) {
            for (int nz = 0; nz < N; nz++) {
                
                for(int i = 0; i < 3; i++) {
                    beta_out(nx, ny, nz, i) = B(nx, ny, nz, i);
                }

            }
        }
    }
}

void Evolution::F_B(SpatialSlice& slice, SpatialSlice& sliceOut) {
    // \partial_t B^i = \partial_t \GammaTilde^i - \eta B^i

    VectorField& B_out = sliceOut.B;

    VectorField& dGammaTilde_dt = sliceOut.GammaTilde; // F(GammaTilde) = dGammaTilde_dt

    VectorField& B = slice.B;

    for (int nx = 0; nx < N; nx++) {
        for (int ny = 0; ny < N; ny++) {
            for (int nz = 0; nz < N; nz++) {
                
                for(int i = 0; i < 3; i++) {
                    B_out(nx, ny, nz, i) = dGammaTilde_dt(nx, ny, nz, i) - (eta * B(nx, ny, nz, i));
                }

            }
        }
    }
}

void Evolution::F_ATilde(SpatialSlice& slice, SpatialSlice& sliceOut) {
    // \partial_t \tilde{A}_{ij} = exp(-4\phi) ( -\partial_j \partial_i \alpha - \Gamma^k_{ji} \partial_k \alpha + \alpha R_ij )^TF + \alpha( K A tilde_ij - 2 A tilde_ik A tilde_lj gamma^kl )

    TensorField& ATilde_out = sliceOut.ATilde;

    ScalarField& phi = slice.phi;
    ScalarField& alpha = slice.alpha;
    TensorField& Ricci = slice.Ricci;
    Tensor3Field& Christoffel = slice.Christoffel;
    TensorField& aux = slice.aux; // to calculate TF
    TensorField& gammaTilde = slice.gammaTilde;
    TensorField& gammaTildeInv = slice.gammaTildeInv;
    ScalarField& K = slice.K;
    TensorField& ATilde = slice.ATilde;

    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                
                //first, aux is calculated
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        aux(nx, ny, nz, i, j) = 0.0;
                        for(int k=0; k<3; k++) {
                            aux(nx, ny, nz, i, j) += -SECOND_PARTIAL_SCALAR(alpha, i, j) - Christoffel(nx, ny, nz, k, j, i) * PARTIAL_SCALAR(alpha, k) + alpha(nx, ny, nz) * Ricci(nx, ny, nz, i, j);
                        }
                    }
                }

                //Now, T_ij ^TF = T_ij - 1/3 gamma_ij gamma^kl T_kl
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        ATilde_out(nx, ny, nz, i, j) = aux(nx, ny, nz, i, j);

                        for(int k = 0; k < 3; k++) {
                            for(int l = 0; l < 3; l++) {
                                ATilde_out(nx, ny, nz, i, j) -= 1.0/3.0 * gammaTilde(nx, ny, nz, i, j) * gammaTildeInv(nx, ny, nz, k, l) * aux(nx, ny, nz, k, l);
                            }
                        }

                        ATilde_out(nx, ny, nz, i, j) *= exp(-4.0 * phi(nx, ny, nz));
                    }
                }

                //Second term: \alpha( K A tilde_ij - 2 A tilde_ik A tilde_lj gamma^kl )
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        ATilde_out(nx, ny, nz, i, j) += alpha(nx, ny, nz) * (K(nx, ny, nz) * ATilde(nx, ny, nz, i, j));

                        for(int k = 0; k < 3; k++) {
                            for(int l = 0; l < 3; l++) {
                                ATilde_out(nx, ny, nz, i, j) -= 2.0 * alpha(nx, ny, nz) * ATilde(nx, ny, nz, i, k) * ATilde(nx, ny, nz, l, j) * gammaTildeInv(nx, ny, nz, k, l);
                            }
                        }
                    }
                }
            }
        }
    }
}

void Evolution::F_GammaTilde(SpatialSlice& slice, SpatialSlice& sliceOut) {
    // \partial_t \tilde{\Gamma}^i = - 2 A^ij gammaInv^ij \partial_j alpha + gammaInv^lj \partial_j \partial_l \beta^i +2alpha( GammaTilde^i_jk ATildeInv^kj -2/3 gammaTildeInv^ij \partial_j K + 6ATildeInv^ij \partial_j phi)

    VectorField& GammaTilde_out = sliceOut.GammaTilde;

    ScalarField& alpha = slice.alpha;
    TensorField& ATildeInv = slice.ATildeInv;
    ScalarField& K = slice.K;
    ScalarField& phi = slice.phi;
    TensorField& gammaInv = slice.gammaInv;
    TensorField& gammaTildeInv = slice.gammaTildeInv;

    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                
                //First term: - 2 A^ij gammaInv^ij \partial_j alpha
                for(int i = 0; i < 3; i++) {
                    GammaTilde_out(nx, ny, nz, i) = 0.0;

                    for(int j = 0; j < 3; j++) {
                        GammaTilde_out(nx, ny, nz, i) -= 2.0 * ATildeInv(nx, ny, nz, i, j) * gammaInv(nx, ny, nz, i, j) * PARTIAL_SCALAR(alpha, j);
                    }
                }

                // Second term: gammaInv^lj \partial_j \partial_l \beta^i
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        for(int l = 0; l < 3; l++) {
                            GammaTilde_out(nx, ny, nz, i) += gammaInv(nx, ny, nz, l, j) * SECOND_PARTIAL_VECTOR(slice.beta, i, j, l);
                        }
                    }
                }

                // Third term: 2alpha( GammaTilde^i_jk ATildeInv^kj )
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        for(int k = 0; k < 3; k++) {
                            GammaTilde_out(nx, ny, nz, i) += 2.0 * alpha(nx, ny, nz) * GammaTilde_out(nx, ny, nz, j) * ATildeInv(nx, ny, nz, k, j);
                        }
                    }
                }

                // Fourth term: -2/3 alpha gammaTildeInv^ij \partial_j K
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        GammaTilde_out(nx, ny, nz, i) -= 2.0/3.0 * gammaTildeInv(nx, ny, nz, i, j) * PARTIAL_SCALAR(K, j) * alpha(nx, ny, nz);
                    }
                }

                // Fifth term: alpha * 6ATildeInv^ij \partial_j phi
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        GammaTilde_out(nx, ny, nz, i) += 6.0 * ATildeInv(nx, ny, nz, i, j) * PARTIAL_SCALAR(phi, j) * alpha(nx, ny, nz);
                    }
                }

            }
        }
    }
}