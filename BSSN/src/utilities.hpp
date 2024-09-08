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

// Define a vector field partial derivative. \partial_i beta^j
#define PARTIAL_VECTOR(beta, i, j) ( \
    (i == 0) ? ( \
        (-beta(PERIODIC_INDEX(nx + 2, NX), ny, nz, j) + 8 * beta(PERIODIC_INDEX(nx + 1, NX), ny, nz, j) - 8 * beta(PERIODIC_INDEX(nx - 1, NX), ny, nz, j) + beta(PERIODIC_INDEX(nx - 2, NX), ny, nz, j)) / (12 * DX)) : \
    (i == 1) ? ( \
        (-beta(nx, PERIODIC_INDEX(ny + 2, NY), nz, j) + 8 * beta(nx, PERIODIC_INDEX(ny + 1, NY), nz, j) - 8 * beta(nx, PERIODIC_INDEX(ny - 1, NY), nz, j) + beta(nx, PERIODIC_INDEX(ny - 2, NY), nz, j)) / (12 * DY)) : \
    (i == 2) ? ( \
        (-beta(nx, ny, PERIODIC_INDEX(nz + 2, NZ), j) + 8 * beta(nx, ny, PERIODIC_INDEX(nz + 1, NZ), j) - 8 * beta(nx, ny, PERIODIC_INDEX(nz - 1, NZ), j) + beta(nx, ny, PERIODIC_INDEX(nz - 2, NZ), j)) / (12 * DZ)) : \
    0.0 \
)

// Define a 2-tensor field partial derivative \partial_i T_jk
#define PARTIAL_TENSOR(T, i, j, k) ( \
    (i == 0) ? ( \
        (-T(PERIODIC_INDEX(nx + 2, NX), ny, nz, j, k) + 8 * T(PERIODIC_INDEX(nx + 1, NX), ny, nz, j, k) - 8 * T(PERIODIC_INDEX(nx - 1, NX), ny, nz, j, k) + T(PERIODIC_INDEX(nx - 2, NX), ny, nz, j, k)) / (12 * DX)) : \
    (i == 1) ? ( \
        (-T(nx, PERIODIC_INDEX(ny + 2, NY), nz, j, k) + 8 * T(nx, PERIODIC_INDEX(ny + 1, NY), nz, j, k) - 8 * T(nx, PERIODIC_INDEX(ny - 1, NY), nz, j, k) + T(nx, PERIODIC_INDEX(ny - 2, NY), nz, j, k)) / (12 * DY)) : \
    (i == 2) ? ( \
        (-T(nx, ny, PERIODIC_INDEX(nz + 2, NZ), j, k) + 8 * T(nx, ny, PERIODIC_INDEX(nz + 1, NZ), j, k) - 8 * T(nx, ny, PERIODIC_INDEX(nz - 1, NZ), j, k) + T(nx, ny, PERIODIC_INDEX(nz - 2, NZ), j, k)) / (12 * DZ)) : \
    0.0 \
)

// \partial_i T_jkl
#define PARTIAL_3TENSOR(T, i, j, k, l) ( \
    (i == 0) ? ( \
        (-T(PERIODIC_INDEX(nx + 2, NX), ny, nz, j, k, l) + 8 * T(PERIODIC_INDEX(nx + 1, NX), ny, nz, j, k, l) - 8 * T(PERIODIC_INDEX(nx - 1, NX), ny, nz, j, k, l) + T(PERIODIC_INDEX(nx - 2, NX), ny, nz, j, k, l)) / (12 * DX)) : \
    (i == 1) ? ( \
        (-T(nx, PERIODIC_INDEX(ny + 2, NY), nz, j, k, l) + 8 * T(nx, PERIODIC_INDEX(ny + 1, NY), nz, j, k, l) - 8 * T(nx, PERIODIC_INDEX(ny - 1, NY), nz, j, k, l) + T(nx, PERIODIC_INDEX(ny - 2, NY), nz, j, k, l)) / (12 * DY)) : \
    (i == 2) ? ( \
        (-T(nx, ny, PERIODIC_INDEX(nz + 2, NZ), j, k, l) + 8 * T(nx, ny, PERIODIC_INDEX(nz + 1, NZ), j, k, l) - 8 * T(nx, ny, PERIODIC_INDEX(nz - 1, NZ), j, k, l) + T(nx, ny, PERIODIC_INDEX(nz - 2, NZ), j, k, l)) / (12 * DZ)) : \
    0.0 \
)

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

//\partial_i \partial_j beta^k
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

// \partial_i \partial_j T_kl
#define SECOND_PARTIAL_TENSOR(T, i, j, k, l) ( \
    (i == 0 && j == 0) ? ( \
        (-T(PERIODIC_INDEX(nx + 2, NX), ny, nz, k, l) + 16 * T(PERIODIC_INDEX(nx + 1, NX), ny, nz, k, l) \
         - 30 * T(nx, ny, nz, k, l) + 16 * T(PERIODIC_INDEX(nx - 1, NX), ny, nz, k, l) - T(PERIODIC_INDEX(nx - 2, NX), ny, nz, k, l)) / (12 * DX * DX)) : \
    (i == 1 && j == 1) ? ( \
        (-T(nx, PERIODIC_INDEX(ny + 2, NY), nz, k, l) + 16 * T(nx, PERIODIC_INDEX(ny + 1, NY), nz, k, l) \
         - 30 * T(nx, ny, nz, k, l) + 16 * T(nx, PERIODIC_INDEX(ny - 1, NY), nz, k, l) - T(nx, PERIODIC_INDEX(ny - 2, NY), nz, k, l)) / (12 * DY * DY)) : \
    (i == 2 && j == 2) ? ( \
        (-T(nx, ny, PERIODIC_INDEX(nz + 2, NZ), k, l) + 16 * T(nx, ny, PERIODIC_INDEX(nz + 1, NZ), k, l) \
         - 30 * T(nx, ny, nz, k, l) + 16 * T(nx, ny, PERIODIC_INDEX(nz - 1, NZ), k, l) - T(nx, ny, PERIODIC_INDEX(nz - 2, NZ), k, l)) / (12 * DZ * DZ)) : \
    (i == 0 && j == 1) ? ( \
        (T(PERIODIC_INDEX(nx + 1, NX), PERIODIC_INDEX(ny + 1, NY), nz, k, l) \
          - T(PERIODIC_INDEX(nx + 1, NX), PERIODIC_INDEX(ny - 1, NY), nz, k, l) \
          - T(PERIODIC_INDEX(nx - 1, NX), PERIODIC_INDEX(ny + 1, NY), nz, k, l) \
            + T(PERIODIC_INDEX(nx - 1, NX), PERIODIC_INDEX(ny - 1, NY), nz, k, l)) / (4 * DX * DY)) : \
    (i == 0 && j == 2) ? ( \
        (T(PERIODIC_INDEX(nx + 1, NX), ny, PERIODIC_INDEX(nz + 1, NZ), k, l) \
          - T(PERIODIC_INDEX(nx + 1, NX), ny, PERIODIC_INDEX(nz - 1, NZ), k, l) \
          - T(PERIODIC_INDEX(nx - 1, NX), ny, PERIODIC_INDEX(nz + 1, NZ), k, l) \
          + T(PERIODIC_INDEX(nx - 1, NX), ny, PERIODIC_INDEX(nz - 1, NZ), k, l)) / (4 * DX * DZ)) : \
    (i == 1 && j == 2) ? ( \
        (T(nx, PERIODIC_INDEX(ny + 1, NY), PERIODIC_INDEX(nz + 1, NZ), k, l) \
          - T(nx, PERIODIC_INDEX(ny + 1, NY), PERIODIC_INDEX(nz - 1, NZ), k, l) \
          - T(nx, PERIODIC_INDEX(ny - 1, NY), PERIODIC_INDEX(nz + 1, NZ), k, l) \
          + T(nx, PERIODIC_INDEX(ny - 1, NY), PERIODIC_INDEX(nz - 1, NZ), k, l)) / (4 * DY * DZ)) : \
    0.0 \
)

// D_i D_j f = partial_i partial_j f - \Gamma^k_ij partial_k f
#define SECOND_COVDERIV_SCALAR(f, i, j, christoffel) ( \
    SECOND_PARTIAL_SCALAR(f, i, j) - christoffel(nx, ny, nz, 0, j, i) * PARTIAL_SCALAR(f, 0) - christoffel(nx, ny, nz, 1, j, i) * PARTIAL_SCALAR(f, 1) - christoffel(nx, ny, nz, 2, j, i) * PARTIAL_SCALAR(f, 2) \
)



#endif