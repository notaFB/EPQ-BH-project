#ifndef FIELD_HPP
#define FIELD_HPP

#include <vector>
#include <array>

#define N 80 // Grid size in each spatial dimension
#define NX N
#define NY N
#define NZ N

#define DX 0.1 // Grid spacing in each spatial dimension
#define DY 0.1
#define DZ 0.1

// Scalar field class
class ScalarField {
public:
    ScalarField() : data(N * N * N, 0.0) {}

    double& operator()(int nx, int ny, int nz) {
        return data[nx * N * N + ny * N + nz];
    }

    const double& operator()(int nx, int ny, int nz) const {
        return data[nx * N * N + ny * N + nz];
    }

private:
    std::vector<double> data;
};


// Vector field class (3 components)
class VectorField {
public:
    VectorField() : data(N * N * N * 3, 0.0) {}

    double& operator()(int nx, int ny, int nz, int i) {
        return data[(nx * N * N + ny * N + nz) * 3 + i];
    }

    const double& operator()(int nx, int ny, int nz, int i) const {
        return data[(nx * N * N + ny * N + nz) * 3 + i];
    }

private:
    std::vector<double> data;
};


// Tensor field class (3x3 components)
class TensorField {
public:
    TensorField() : data(N * N * N * 9, 0.0) {}

    double& operator()(int nx, int ny, int nz, int i, int j) {
        return data[(nx * N * N + ny * N + nz) * 9 + i * 3 + j];
    }

    const double& operator()(int nx, int ny, int nz, int i, int j) const {
        return data[(nx * N * N + ny * N + nz) * 9 + i * 3 + j];
    }

private:
    std::vector<double> data;
};


// Tensor field class with 3 indices (3x3x3 components)
class Tensor3Field {
public:
    Tensor3Field() : data(N * N * N * 27, 0.0) {}

    double& operator()(int nx, int ny, int nz, int i, int j, int k) {
        return data[(nx * N * N + ny * N + nz) * 27 + i * 9 + j * 3 + k];
    }

    const double& operator()(int nx, int ny, int nz, int i, int j, int k) const {
        return data[(nx * N * N + ny * N + nz) * 27 + i * 9 + j * 3 + k];
    }

private:
    std::vector<double> data;
};

// Spatial slice class containing all relevant fields
class SpatialSlice {
public:
    ScalarField phi;
    ScalarField K;

    ScalarField alpha;
    VectorField beta;
    VectorField betaInv;

    TensorField gammaTilde;
    TensorField gammaTildeInv;

    TensorField gamma;
    TensorField gammaInv;
    
    TensorField ATilde;
    TensorField ATildeInv;

    Tensor3Field Christoffel;
    Tensor3Field ChristoffelTilde;
    TensorField Ricci;
};

#endif // FIELD_HPP
