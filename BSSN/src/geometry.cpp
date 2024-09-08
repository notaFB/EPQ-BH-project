
#include "geometry.hpp"
#include <cmath>

void Geometry::computegamma(SpatialSlice& slice)
{
    //gamma_ij = exp(4 phi) gammaTilde_ij
    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        slice.gamma(nx, ny, nz, i, j) = exp(4.0 * slice.phi(nx, ny, nz)) * slice.gammaTilde(nx, ny, nz, i, j);
                    }
                }
            }
        }
    }
}

void Geometry::computeInverseMetric(TensorField& gamma, TensorField& gammaInv) {
    // Compute the inverse metric

    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                
                //Compute the determinant of the metric
                double det = gamma(nx, ny, nz, 0, 0) * (gamma(nx, ny, nz, 1, 1) * gamma(nx, ny, nz, 2, 2) - gamma(nx, ny, nz, 1, 2) * gamma(nx, ny, nz, 2, 1)) - gamma(nx, ny, nz, 0, 1) * (gamma(nx, ny, nz, 1, 0) * gamma(nx, ny, nz, 2, 2) - gamma(nx, ny, nz, 1, 2) * gamma(nx, ny, nz, 2, 0)) + gamma(nx, ny, nz, 0, 2) * (gamma(nx, ny, nz, 1, 0) * gamma(nx, ny, nz, 2, 1) - gamma(nx, ny, nz, 1, 1) * gamma(nx, ny, nz, 2, 0));

                //Compute the inverse metric
                gammaInv(nx, ny, nz, 0, 0) = (gamma(nx, ny, nz, 1, 1) * gamma(nx, ny, nz, 2, 2) - gamma(nx, ny, nz, 1, 2) * gamma(nx, ny, nz, 2, 1)) / det;
                gammaInv(nx, ny, nz, 0, 1) = (gamma(nx, ny, nz, 0, 2) * gamma(nx, ny, nz, 2, 1) - gamma(nx, ny, nz, 0, 1) * gamma(nx, ny, nz, 2, 2)) / det;
                gammaInv(nx, ny, nz, 0, 2) = (gamma(nx, ny, nz, 0, 1) * gamma(nx, ny, nz, 1, 2) - gamma(nx, ny, nz, 0, 2) * gamma(nx, ny, nz, 1, 1)) / det;
                gammaInv(nx, ny, nz, 1, 0) = (gamma(nx, ny, nz, 1, 2) * gamma(nx, ny, nz, 2, 0) - gamma(nx, ny, nz, 1, 0) * gamma(nx, ny, nz, 2, 2)) / det;
                gammaInv(nx, ny, nz, 1, 1) = (gamma(nx, ny, nz, 0, 0) * gamma(nx, ny, nz, 2, 2) - gamma(nx, ny, nz, 0, 2) * gamma(nx, ny, nz, 2, 0)) / det;
                gammaInv(nx, ny, nz, 1, 2) = (gamma(nx, ny, nz, 0, 2) * gamma(nx, ny, nz, 1, 0) - gamma(nx, ny, nz, 0, 0) * gamma(nx, ny, nz, 1, 2)) / det;
                gammaInv(nx, ny, nz, 2, 0) = (gamma(nx, ny, nz, 1, 0) * gamma(nx, ny, nz, 2, 1) - gamma(nx, ny, nz, 1, 1) * gamma(nx, ny, nz, 2, 0)) / det;
                gammaInv(nx, ny, nz, 2, 1) = (gamma(nx, ny, nz, 0, 1) * gamma(nx, ny, nz, 2, 0) - gamma(nx, ny, nz, 0, 0) * gamma(nx, ny, nz, 2, 1)) / det;
                gammaInv(nx, ny, nz, 2, 2) = (gamma(nx, ny, nz, 0, 0) * gamma(nx, ny, nz, 1, 1) - gamma(nx, ny, nz, 0, 1) * gamma(nx, ny, nz, 1, 0)) / det;

            }
        }
    }
}

void Geometry::computeInverseTensor(TensorField& T, TensorField& TInv, TensorField& gammaInv) {
    // Compute the inverse of a tensor
    //T^ij = gamma^ik gamma^jl T_kl

    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        TInv(nx, ny, nz, i, j) = 0.0;
                        for(int k = 0; k < 3; k++) {
                            for(int l = 0; l < 3; l++) {
                                TInv(nx, ny, nz, i, j) += gammaInv(nx, ny, nz, i, k) * gammaInv(nx, ny, nz, j, l) * T(nx, ny, nz, k, l);
                            }
                        }
                    }
                }

            }
        }
    }

}

void Geometry::computeChristoffel(TensorField& gamma, TensorField& gammaInv, Tensor3Field& Christoffel) {
    // Compute the Christoffel symbols of the second kind
    //Gamma^i_jk = 1/2 gamma^il (partial_j gamma_lk + partial_k gamma_jl - partial_l gamma_jk)

    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        for(int k = 0; k < 3; k++) {

                            Christoffel(nx, ny, nz, i, j, k) = 0.0;
                            for(int l = 0; l < 3; l++) {
                                Christoffel(nx, ny, nz, i, j, k) += 0.5 * gammaInv(nx, ny, nz, i, l) * (PARTIAL_TENSOR(gamma, j, l, k) + PARTIAL_TENSOR(gamma, k, l, j) - PARTIAL_TENSOR(gamma, l, j, k));
                            }
                        }
                    }
                }

            }
        }
    }

}

void Geometry::computerRicci(SpatialSlice& slice) {
    // Compute the Ricci tensor

    
    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        slice.Ricci(nx, ny, nz, i, j) = 0.0;

                        //First term: -2 D tilde_i D tilde_j phi 
                        slice.Ricci(nx, ny, nz, i, j) -= 2.0 * SECOND_COVDERIV_SCALAR(slice.phi, i, j, slice.ChristoffelTilde);

                        //Second term: -2 gamma tilde_ij gamma tilde^lk D_k D_l phi
                        for(int k = 0; k < 3; k++) {
                            for(int l = 0; l < 3; l++) {
                                slice.Ricci(nx, ny, nz, i, j) -= 2.0 * slice.gammaTilde(nx, ny, nz, i, j) * slice.gammaTildeInv(nx, ny, nz, k, l) * SECOND_COVDERIV_SCALAR(slice.phi, k, l, slice.ChristoffelTilde);
                            }
                        }

                        //Third term: +4 partial_i phi partial_j phi
                        slice.Ricci(nx, ny, nz, i, j) += 4.0 * PARTIAL_SCALAR(slice.phi, i) * PARTIAL_SCALAR(slice.phi, j);

                        //Fourth term: -4 gamma tilde_ij gamma tilde^kl partial_k phi partial_l phi
                        for(int k = 0; k < 3; k++) {
                            for(int l = 0; l < 3; l++) {
                                slice.Ricci(nx, ny, nz, i, j) -= 4.0 * slice.gammaTilde(nx, ny, nz, i, j) * slice.gammaTildeInv(nx, ny, nz, k, l) * PARTIAL_SCALAR(slice.phi, k) * PARTIAL_SCALAR(slice.phi, l);
                            }
                        }


                        //Now second half: tilde R_ij

                        //First term: -1/2 gamma tilde^lm \partial_m \partial_l gamma tilde_ij
                        for(int l = 0; l < 3; l++) {
                            for(int m = 0; m < 3; m++) {
                                slice.Ricci(nx, ny, nz, i, j) -= 0.5 * slice.gammaTildeInv(nx, ny, nz, l, m) * SECOND_PARTIAL_TENSOR(slice.gammaTilde, m, l, i, j);
                            }
                        }

                        //Seecond term: +\tilde \gamma_{ki} \partial_j \tilde \Gamma^k +\tilde \gamma_{kj} \partial_i \tilde \Gamma^k
                        for(int k = 0; k < 3; k++) {
                            slice.Ricci(nx, ny, nz, i, j) += slice.gammaTilde(nx, ny, nz, k, i) * PARTIAL_VECTOR(slice.GammaTilde, j, k) + slice.gammaTilde(nx, ny, nz, k, j) * PARTIAL_VECTOR(slice.GammaTilde, i, k);
                        }

                        //Third term: +\tilde \Gamma^k ( \tilde \gamma_{li} \tilde \Gamma^l_{jk} + \tilde \gamma_{lj} \tilde \Gamma^l_{ik} )
                        for(int k = 0; k < 3; k++) {
                            for(int l = 0; l < 3; l++) {
                                slice.Ricci(nx, ny, nz, i, j) += slice.GammaTilde(nx, ny, nz, k) * (slice.gammaTilde(nx, ny, nz, l, i) * slice.ChristoffelTilde(nx, ny, nz, l, j, k) + slice.gammaTilde(nx, ny, nz, l, j) * slice.ChristoffelTilde(nx, ny, nz, l, i, k));
                            }
                        }

                        //Fourth term: +\tilde \gamma^{lm} ( 2\tilde \Gamma^k_{li} \tilde \gamma_{nj} \tilde \Gamma^n_{km} + 2\tilde \Gamma^k_{lj} \tilde \gamma_{ni} \tilde \Gamma^n{km} )
                        for(int l = 0; l < 3; l++) {
                            for(int m = 0; m < 3; m++) {
                                for(int n = 0; n < 3; n++) {
                                    for(int k = 0; k < 3; k++) {
                                        slice.Ricci(nx, ny, nz, i, j) += slice.gammaTildeInv(nx, ny, nz, l, m) * (2.0 * slice.ChristoffelTilde(nx, ny, nz, k, l, i) * slice.gammaTilde(nx, ny, nz, n, j) * slice.ChristoffelTilde(nx, ny, nz, n, k, m) + 2.0 * slice.ChristoffelTilde(nx, ny, nz, k, l, j) * slice.gammaTilde(nx, ny, nz, n, i) * slice.ChristoffelTilde(nx, ny, nz, n, k, m));
                                    }
                                }
                            }
                        }

                        //Last term: + \tilde \gamma^{lm} \tilde \gamma_{nk} \tilde \Gamma^k_{im} \tilde \Gamma^n_{lj}
                        for(int l = 0; l < 3; l++) {
                            for(int m = 0; m < 3; m++) {
                                for(int n = 0; n < 3; n++) {
                                    for(int k = 0; k < 3; k++) {
                                        slice.Ricci(nx, ny, nz, i, j) += slice.gammaTildeInv(nx, ny, nz, l, m) * slice.gammaTilde(nx, ny, nz, n, k) * slice.ChristoffelTilde(nx, ny, nz, k, i, m) * slice.ChristoffelTilde(nx, ny, nz, n, j, l);
                                    }
                                }
                            }
                        }

                    }
                }

            }
        }
    }
}

void Geometry::computeGeometry(SpatialSlice& slice) {
    // Compute all the geometry

    computegamma(slice);

    computeInverseMetric(slice.gamma, slice.gammaInv);
    computeInverseMetric(slice.gammaTilde, slice.gammaTildeInv);
    computeInverseTensor(slice.ATilde, slice.ATildeInv, slice.gammaTildeInv);
    computeChristoffel(slice.gamma, slice.gammaInv, slice.Christoffel);
    computeChristoffel(slice.gammaTilde, slice.gammaTildeInv, slice.ChristoffelTilde);

    computerRicci(slice);

}