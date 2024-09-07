#ifndef UTILITIES_HPP
#define UTILITIES_HPP

// Macro for periodic boundary handling
#define PERIODIC_INDEX(i, N) (((i) + (N)) % (N))

// Define the partial derivative macro with periodic boundary conditions
#define PARTIAL_SCALAR(alpha, i) ( \
    (i == 0) ? (-alpha(PERIODIC_INDEX(nx + 2, NX), ny, nz) + 8 * alpha(PERIODIC_INDEX(nx + 1, NX), ny, nz) - 8 * alpha(PERIODIC_INDEX(nx - 1, NX), ny, nz) + alpha(PERIODIC_INDEX(nx - 2, NX), ny, nz)) / (12 * DX) : \
    (i == 1) ? (-alpha(nx, PERIODIC_INDEX(ny + 2, NY), nz) + 8 * alpha(nx, PERIODIC_INDEX(ny + 1, NY), nz) - 8 * alpha(nx, PERIODIC_INDEX(ny - 1, NY), nz) + alpha(nx, PERIODIC_INDEX(ny - 2, NY), nz)) / (12 * DY) : \
    (i == 2) ? (-alpha(nx, ny, PERIODIC_INDEX(nz + 2, NZ)) + 8 * alpha(nx, ny, PERIODIC_INDEX(nz + 1, NZ)) - 8 * alpha(nx, ny, PERIODIC_INDEX(nz - 1, NZ)) + alpha(nx, ny, PERIODIC_INDEX(nz - 2, NZ))) / (12 * DZ) : \
    0.0)

// Define a vector field partial derivative
#define PARTIAL_VECTOR(alpha, i, j) ( \
    (i == 0) ? ( \
        (j == 0) ? (-alpha(PERIODIC_INDEX(nx + 2, NX), ny, nz) + 8 * alpha(PERIODIC_INDEX(nx + 1, NX), ny, nz) - 8 * alpha(PERIODIC_INDEX(nx - 1, NX), ny, nz) + alpha(PERIODIC_INDEX(nx - 2, NX), ny, nz)) / (12 * DX) : \
        (j == 1) ? (-alpha(nx, PERIODIC_INDEX(ny + 2, NY), nz) + 8 * alpha(nx, PERIODIC_INDEX(ny + 1, NY), nz) - 8 * alpha(nx, PERIODIC_INDEX(ny - 1, NY), nz) + alpha(nx, PERIODIC_INDEX(ny - 2, NY), nz)) / (12 * DY) : \
        (j == 2) ? (-alpha(nx, ny, PERIODIC_INDEX(nz + 2, NZ)) + 8 * alpha(nx, ny, PERIODIC_INDEX(nz + 1, NZ)) - 8 * alpha(nx, ny, PERIODIC_INDEX(nz - 1, NZ)) + alpha(nx, ny, PERIODIC_INDEX(nz - 2, NZ))) / (12 * DZ) : \
        0.0 \
    ) : \
    (i == 1) ? ( \
        (j == 0) ? (-alpha(PERIODIC_INDEX(nx + 2, NX), ny, nz) + 8 * alpha(PERIODIC_INDEX(nx + 1, NX), ny, nz) - 8 * alpha(PERIODIC_INDEX(nx - 1, NX), ny, nz) + alpha(PERIODIC_INDEX(nx - 2, NX), ny, nz)) / (12 * DX) : \
        (j == 1) ? (-alpha(nx, PERIODIC_INDEX(ny + 2, NY), nz) + 8 * alpha(nx, PERIODIC_INDEX(ny + 1, NY), nz) - 8 * alpha(nx, PERIODIC_INDEX(ny - 1, NY), nz) + alpha(nx, PERIODIC_INDEX(ny - 2, NY), nz)) / (12 * DY) : \
        (j == 2) ? (-alpha(nx, ny, PERIODIC_INDEX(nz + 2, NZ)) + 8 * alpha(nx, ny, PERIODIC_INDEX(nz + 1, NZ)) - 8 * alpha(nx, ny, PERIODIC_INDEX(nz - 1, NZ)) + alpha(nx, ny, PERIODIC_INDEX(nz - 2, NZ))) / (12 * DZ) : \
        0.0 \
    ) : \
    0.0)

// Define second-order mixed partial derivative macro for scalar field with periodic boundary conditions
#define SECOND_PARTIAL_SCALAR(alpha, i, j) ( \
    (i == 0 && j == 0) ? ( \
        (-alpha(PERIODIC_INDEX(nx + 2, NX), ny, nz) + 16 * alpha(PERIODIC_INDEX(nx + 1, NX), ny, nz) \
         - 30 * alpha(nx, ny, nz) + 16 * alpha(PERIODIC_INDEX(nx - 1, NX), ny, nz) - alpha(PERIODIC_INDEX(nx - 2, NX), ny, nz)) / (12 * DX * DX)) : \
    (i == 1 && j == 1) ? ( \
        (-alpha(nx, PERIODIC_INDEX(ny + 2, NY), nz) + 16 * alpha(nx, PERIODIC_INDEX(ny + 1, NY), nz) \
         - 30 * alpha(nx, ny, nz) + 16 * alpha(nx, PERIODIC_INDEX(ny - 1, NY), nz) - alpha(nx, PERIODIC_INDEX(ny - 2, NY), nz)) / (12 * DY * DY)) : \
    (i == 2 && j == 2) ? ( \
        (-alpha(nx, ny, PERIODIC_INDEX(nz + 2, NZ)) + 16 * alpha(nx, ny, PERIODIC_INDEX(nz + 1, NZ)) \
         - 30 * alpha(nx, ny, nz) + 16 * alpha(nx, ny, PERIODIC_INDEX(nz - 1, NZ)) - alpha(nx, ny, PERIODIC_INDEX(nz - 2, NZ))) / (12 * DZ * DZ)) : \
    (i == 0 && j == 1) ? ( \
        (-(alpha(PERIODIC_INDEX(nx + 1, NX), PERIODIC_INDEX(ny + 1, NY), nz) \
          - alpha(PERIODIC_INDEX(nx + 1, NX), PERIODIC_INDEX(ny - 1, NY), nz)) \
          - (alpha(PERIODIC_INDEX(nx - 1, NX), PERIODIC_INDEX(ny + 1, NY), nz) \
          - alpha(PERIODIC_INDEX(nx - 1, NX), PERIODIC_INDEX(ny - 1, NY), nz))) / (4 * DX * DY)) : \
    (i == 0 && j == 2) ? ( \
        (-(alpha(PERIODIC_INDEX(nx + 1, NX), ny, PERIODIC_INDEX(nz + 1, NZ)) \
          - alpha(PERIODIC_INDEX(nx + 1, NX), ny, PERIODIC_INDEX(nz - 1, NZ))) \
          - (alpha(PERIODIC_INDEX(nx - 1, NX), ny, PERIODIC_INDEX(nz + 1, NZ)) \
          - alpha(PERIODIC_INDEX(nx - 1, NX), ny, PERIODIC_INDEX(nz - 1, NZ)))) / (4 * DX * DZ)) : \
    (i == 1 && j == 2) ? ( \
        (-(alpha(nx, PERIODIC_INDEX(ny + 1, NY), PERIODIC_INDEX(nz + 1, NZ)) \
          - alpha(nx, PERIODIC_INDEX(ny + 1, NY), PERIODIC_INDEX(nz - 1, NZ))) \
          - (alpha(nx, PERIODIC_INDEX(ny - 1, NY), PERIODIC_INDEX(nz + 1, NZ)) \
          - alpha(nx, PERIODIC_INDEX(ny - 1, NY), PERIODIC_INDEX(nz - 1, NZ)))) / (4 * DY * DZ)) : \
    0.0 \
)

#define SECOND_PARTIAL_VECTOR(beta, k, i, j) ( \
    (i == 0 && j == 0) ? ( \
        (-beta(PERIODIC_INDEX(nx + 2, NX), ny, nz, k) + 16 * beta(PERIODIC_INDEX(nx + 1, NX), ny, nz, k) \
         - 30 * beta(nx, ny, nz, k) + 16 * beta(PERIODIC_INDEX(nx - 1, NX), ny, nz, k) - beta(PERIODIC_INDEX(nx - 2, NX), ny, nz, k)) / (12 * DX * DX)) : \
    (i == 1 && j == 1) ? ( \
        (-beta(nx, PERIODIC_INDEX(ny + 2, NY), nz, k) + 16 * beta(nx, PERIODIC_INDEX(ny + 1, NY), nz, k) \
         - 30 * beta(nx, ny, nz, k) + 16 * beta(nx, PERIODIC_INDEX(ny - 1, NY), nz, k) - beta(nx, PERIODIC_INDEX(ny - 2, NY), nz, k)) / (12 * DY * DY)) : \
    (i == 2 && j == 2) ? ( \
        (-beta(nx, ny, PERIODIC_INDEX(nz + 2, NZ), k) + 16 * beta(nx, ny, PERIODIC_INDEX(nz + 1, NZ), k) \
         - 30 * beta(nx, ny, nz, k) + 16 * beta(nx, ny, PERIODIC_INDEX(nz - 1, NZ), k) - beta(nx, ny, PERIODIC_INDEX(nz - 2, NZ), k)) / (12 * DZ * DZ)) : \
    (i == 0 && j == 1) ? ( \
        (beta(PERIODIC_INDEX(nx + 1, NX), PERIODIC_INDEX(ny + 1, NY), nz, k) \
          - beta(PERIODIC_INDEX(nx + 1, NX), PERIODIC_INDEX(ny - 1, NY), nz, k) \
          - beta(PERIODIC_INDEX(nx - 1, NX), PERIODIC_INDEX(ny + 1, NY), nz, k) \
          + beta(PERIODIC_INDEX(nx - 1, NX), PERIODIC_INDEX(ny - 1, NY), nz, k)) / (4 * DX * DY)) : \
    (i == 0 && j == 2) ? ( \
        (beta(PERIODIC_INDEX(nx + 1, NX), ny, PERIODIC_INDEX(nz + 1, NZ), k) \
          - beta(PERIODIC_INDEX(nx + 1, NX), ny, PERIODIC_INDEX(nz - 1, NZ), k) \
          - beta(PERIODIC_INDEX(nx - 1, NX), ny, PERIODIC_INDEX(nz + 1, NZ), k) \
          + beta(PERIODIC_INDEX(nx - 1, NX), ny, PERIODIC_INDEX(nz - 1, NZ), k)) / (4 * DX * DZ)) : \
    (i == 1 && j == 2) ? ( \
        (beta(nx, PERIODIC_INDEX(ny + 1, NY), PERIODIC_INDEX(nz + 1, NZ), k) \
          - beta(nx, PERIODIC_INDEX(ny + 1, NY), PERIODIC_INDEX(nz - 1, NZ), k) \
          - beta(nx, PERIODIC_INDEX(ny - 1, NY), PERIODIC_INDEX(nz + 1, NZ), k) \
          + beta(nx, PERIODIC_INDEX(ny - 1, NY), PERIODIC_INDEX(nz - 1, NZ), k)) / (4 * DY * DZ)) : \
    0.0 \
)

// D_i D_j f = partial_i partial_j f - \Gamma^k_ij partial_k f
#define SECOND_COVDERIV_SCALAR(f, i, j, christoffel) ( \
    SECOND_PARTIAL_SCALAR(f, i, j) - christoffel(nx, ny, nz, 0, j, i) * PARTIAL_SCALAR(f, 0) - christoffel(nx, ny, nz, 1, j, i) * PARTIAL_SCALAR(f, 1) - christoffel(nx, ny, nz, 2, j, i) * PARTIAL_SCALAR(f, 2) \
)



#endif