
#include "geometry.hpp"

void Geometry::computeInverseMetric(TensorField& gamma, TensorField& gammaInv) {
    // Compute the inverse metric
    (void) gamma;
    (void) gammaInv;

}

void Geometry::computeInverseTensor(TensorField& T, TensorField& TInv, TensorField& gammaInv) {
    // Compute the inverse of a tensor
    (void) T;
    (void) TInv;
    (void) gammaInv;

}

void Geometry::computeChristoffel(TensorField& gamma, TensorField& gammaInv, Tensor3Field& Christoffel) {
    // Compute the Christoffel symbols of the second kind
    (void) gamma;
    (void) gammaInv;
    (void) Christoffel;

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

                    }
                }

            }
        }
    }

}

void Geometry::compute(SpatialSlice& slice) {
    // Compute all the geometry

    computeInverseMetric(slice.gamma, slice.gammaInv);
    computeInverseMetric(slice.gammaTilde, slice.gammaTildeInv);
    computeInverseTensor(slice.ATilde, slice.ATildeInv, slice.gammaTildeInv);
    computeChristoffel(slice.gamma, slice.gammaInv, slice.Christoffel);
    computeChristoffel(slice.gammaTilde, slice.gammaTildeInv, slice.ChristoffelTilde);

    computerRicci(slice);

}