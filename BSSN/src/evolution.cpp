
#include"evolution.hpp"
#include "utilities.hpp"

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
    TensorField& gammaTildeInv = slice.gammaTildeInv;

    for (int nx = 0; nx < N; nx++) {
        for (int ny = 0; ny < N; ny++) {
            for (int nz = 0; nz < N; nz++) {
                
                K_out(nx, ny, nz) = 0.0;

                //First term: gamma^{ij} \left( \partial_j \partial_i \alpha - \Gamma^k_{ji} \partial_k \alpha \right)
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {
                        K_out(nx, ny, nz) += gammaInv(nx, ny, nz, i, j) * (SECOND_PARTIAL_SCALAR(alpha, i, j) - Christoffel(nx, ny, nz, j, i, 0) * PARTIAL_SCALAR(alpha, 0) - Christoffel(nx, ny, nz, j, i, 1) * PARTIAL_SCALAR(alpha, 1) - Christoffel(nx, ny, nz, j, i, 2) * PARTIAL_SCALAR(alpha, 2));
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
                
            }
        }
    }
}