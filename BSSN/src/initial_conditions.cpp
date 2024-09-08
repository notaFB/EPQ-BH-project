
#include "initial_conditions.hpp"

#include <cmath>

void InitialConditions::setInitialConditions(SpatialSlice& slice) {
    // Set the initial conditions for the simulation

    setgammaTilde(slice);
    setATilde(slice);
    setK(slice);
    setGammaTilde(slice);
    setPhi(slice);
    setalpha(slice);
    setBeta(slice);
    setB(slice);

}

void InitialConditions::setgammaTilde(SpatialSlice& slice) {
    // Set the initial conformal metric

    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                slice.gammaTilde(nx, ny, nz, 0, 0) = 1.0;
                slice.gammaTilde(nx, ny, nz, 1, 1) = 1.0;
                slice.gammaTilde(nx, ny, nz, 2, 2) = 1.0;
                slice.gammaTilde(nx, ny, nz, 0, 1) = 0.0;
                slice.gammaTilde(nx, ny, nz, 0, 2) = 0.0;
                slice.gammaTilde(nx, ny, nz, 1, 2) = 0.0;
            }
        }
    }

}

void InitialConditions::setATilde(SpatialSlice& slice) {
    // Set the initial conformal extrinsic curvature

    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                slice.ATilde(nx, ny, nz, 0, 0) = 0.0;
                slice.ATilde(nx, ny, nz, 1, 1) = 0.0;
                slice.ATilde(nx, ny, nz, 2, 2) = 0.0;
                slice.ATilde(nx, ny, nz, 0, 1) = 0.0;
                slice.ATilde(nx, ny, nz, 0, 2) = 0.0;
                slice.ATilde(nx, ny, nz, 1, 2) = 0.0;
            }
        }
    }

}

void InitialConditions::setK(SpatialSlice& slice) {
    // Set the initial trace of the extrinsic curvature

    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                slice.K(nx, ny, nz) = 0.0;
            }
        }
    }

}

void InitialConditions::setGammaTilde(SpatialSlice& slice) {
    // Set to zero

    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                slice.GammaTilde(nx, ny, nz, 0) = 0.0;
                slice.GammaTilde(nx, ny, nz, 1) = 0.0;
                slice.GammaTilde(nx, ny, nz, 2) = 0.0;
            }
        }
    }

}

void InitialConditions::setPhi(SpatialSlice& slice) {
    // phi = ln ( 1 + M/2r ), where r=sqrt((x-x_c)^2 + (y-y_c)^2 + (z-z_c)^2)
    double L = (NX * DX) / 2.0;
    double M = 0.05;
    double x_c = L / 2.0;
    double y_c = L / 2.0;
    double z_c = L / 2.0;

    double p1_x = x_c - (L/8.0);
    double p1_y = y_c - (L/8.0);

    double p2_x = x_c + (L/8.0);
    double p2_y = y_c + (L/8.0);


    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {

                double x1 = (nx * DX)/2.0 - p1_x;
                double y1 = (ny * DY)/2.0 - p1_y;
                double z1 = (nz * DZ)/2.0 - z_c;

                double r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);

                double x2 = (nx * DX)/2.0 - p2_x;
                double y2 = (ny * DY)/2.0 - p2_y;
                double z2 = (nz * DZ)/2.0 - z_c;

                double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);

                slice.phi(nx, ny, nz) = log(1.0 + M / (2.0 * r1) + M / (2.0 * r2));
                
            }
        }
    }

}

void InitialConditions::setalpha(SpatialSlice& slice) {
    // alpha = ( 1 - M/2r ) / (1 + M/2r)
    double L = (NX * DX) / 2.0;
    double M = 0.05;
    double x_c = L / 2.0;
    double y_c = L / 2.0;
    double z_c = L / 2.0;

    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {

                double x = (nx * DX)/2.0 - x_c;
                double y = (ny * DY)/2.0 - y_c;
                double z = (nz * DZ)/2.0 - z_c;

                double r = sqrt(x * x + y * y + z * z);

                slice.alpha(nx, ny, nz) = (1.0 - M / (2.0 * r)) / (1.0 + M / (2.0 * r));
            }
        }
    }

}

void InitialConditions::setBeta(SpatialSlice& slice) {
    // beta = 0
    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                slice.beta(nx, ny, nz, 0) = 0.0;
                slice.beta(nx, ny, nz, 1) = 0.0;
                slice.beta(nx, ny, nz, 2) = 0.0;
            }
        }
    }

}

void InitialConditions::setB(SpatialSlice& slice) {
    // B = 0
    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                slice.B(nx, ny, nz, 0) = 0.0;
                slice.B(nx, ny, nz, 1) = 0.0;
                slice.B(nx, ny, nz, 2) = 0.0;
            }
        }
    }

}