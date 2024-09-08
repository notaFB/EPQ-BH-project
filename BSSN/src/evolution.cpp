
#include"evolution.hpp"
#include "utilities.hpp"
#include "geometry.hpp"
#include <cmath>
#include<iostream>

void Evolution::RK4_step()
{

    

    //Get the data from FieldData
    FieldData& fieldData = FieldData::getInstance();
    SpatialSlice& v1 = fieldData.getV1Slice(); // v1 stores the current state vector s
    SpatialSlice& k1 = fieldData.getk1Slice(); // k1 stores the first derivative of s
    SpatialSlice& k2 = fieldData.getk2Slice(); // k2 stores the second derivative of s
    SpatialSlice& k3 = fieldData.getk3Slice(); // k3 stores the third derivative of s
    SpatialSlice& k4 = fieldData.getk4Slice(); // k4 stores the fourth derivative of s
    SpatialSlice& aux = fieldData.getAuxSlice(); // aux is used as a temporary storage

    double halfStep = 0.5;
    double fullStep = 1.0;

    Geometry::computeGeometry(v1);
    //cout 0,0,0 value of phi
    //std::cout << "phi(0,0,0) = " << v1.phi(0,0,0) << std::endl;
    //std::cout << "alpha(0,0,0) = " << v1.alpha(0,0,0) << std::endl;

    //First step: k1 = F(v1)
    Evolution::F(v1, k1);
    Geometry::computeGeometry(k1);
    //cout 0,0,0 value of phi
    //std::cout << "k1 phi(0,0,0) = " << k1.phi(0,0,0) << std::endl;
    //std::cout << "k1 alpha(0,0,0) = " << k1.alpha(0,0,0) << std::endl;

    //Second step: k2 = F(v1 + 0.5*dt*k1). Use aux as a temporary storage
    aux = k1;
    //std::cout << "aux phi(0,0,0) = " << aux.phi(0,0,0) << std::endl;
    //std::cout << "aux alpha(0,0,0) = " << aux.alpha(0,0,0) << std::endl;
    aux *= halfStep * DT;
    aux += v1;
    Geometry::computeGeometry(aux);
    //std::cout << "aux' phi(0,0,0) = " << aux.phi(0,0,0) << std::endl;
    //std::cout << "aux' alpha(0,0,0) = " << aux.alpha(0,0,0) << std::endl;
    Evolution::F(aux, k2);
    Geometry::computeGeometry(k2);
    //std::cout << "k2 phi(0,0,0) = " << k2.phi(0,0,0) << std::endl;
    //std::cout << "k2 alpha(0,0,0) = " << k2.alpha(0,0,0) << std::endl;

    //Third step: k3 = F(v1 + 0.5*dt*k2)
    aux = k2;
    aux *= halfStep * DT;
    aux += v1;
    Geometry::computeGeometry(aux);
    Evolution::F(aux, k3);
    Geometry::computeGeometry(k3);
    //std::cout << "k3 phi(0,0,0) = " << k3.phi(0,0,0) << std::endl;

    //Fourth step: k4 = F(v1 + dt*k3)
    aux = k3;
    aux *= fullStep * DT;
    aux += v1;
    Geometry::computeGeometry(aux);
    Evolution::F(aux, k4);
    Geometry::computeGeometry(k4);
    //std::cout << "k4 phi(0,0,0) = " << k4.phi(0,0,0) << std::endl;

    //Now update v1. v1 = v1 + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    k1 *= DT / 6.0;
    k2 *= 2.0 * DT / 6.0;
    k3 *= 2.0 * DT / 6.0;
    k4 *= DT / 6.0;

    v1 += k1;
    v1 += k2;
    v1 += k3;
    v1 += k4;

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
                            GammaTilde_out(nx, ny, nz, i) += gammaInv(nx, ny, nz, l, j) * SECOND_PARTIAL_VECTOR(slice.beta, j, l, i);
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