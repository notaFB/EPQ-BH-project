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

#define eta 1.0 // Damping parameter,  ~ M

// Base class template for field operations
template <typename Derived>
class FieldBase {
public:
    Derived& operator+=(const Derived& rhs) {
        for (size_t i = 0; i < data.size(); ++i) {
            data[i] += rhs.data[i];
        }
        return static_cast<Derived&>(*this);
    }

    Derived& operator*=(const double& rhs) {
        for (size_t i = 0; i < data.size(); ++i) {
            data[i] *= rhs;
        }
        return static_cast<Derived&>(*this);
    }

    Derived& operator=(const Derived& rhs) {
        for (size_t i = 0; i < data.size(); ++i) {
            data[i] = rhs.data[i];
        }
        return static_cast<Derived&>(*this);
    }

protected:
    std::vector<double> data;
};

// Scalar field class
class ScalarField : public FieldBase<ScalarField> {
public:
    ScalarField() { data.resize(N * N * N, 0.0); }

    double& operator()(int nx, int ny, int nz) {
        return data[nx * N * N + ny * N + nz];
    }

    const double& operator()(int nx, int ny, int nz) const {
        return data[nx * N * N + ny * N + nz];
    }
};

// Vector field class (3 components)
class VectorField : public FieldBase<VectorField> {
public:
    VectorField() { data.resize(N * N * N * 3, 0.0); }

    double& operator()(int nx, int ny, int nz, int i) {
        return data[(nx * N * N + ny * N + nz) * 3 + i];
    }

    const double& operator()(int nx, int ny, int nz, int i) const {
        return data[(nx * N * N + ny * N + nz) * 3 + i];
    }
};

// Tensor field class (3x3 components)
class TensorField : public FieldBase<TensorField> {
public:
    TensorField() { data.resize(N * N * N * 9, 0.0); }

    double& operator()(int nx, int ny, int nz, int i, int j) {
        return data[(nx * N * N + ny * N + nz) * 9 + i * 3 + j];
    }

    const double& operator()(int nx, int ny, int nz, int i, int j) const {
        return data[(nx * N * N + ny * N + nz) * 9 + i * 3 + j];
    }
};

// Tensor field class with 3 indices (3x3x3 components)
class Tensor3Field : public FieldBase<Tensor3Field> {
public:
    Tensor3Field() { data.resize(N * N * N * 27, 0.0); }

    double& operator()(int nx, int ny, int nz, int i, int j, int k) {
        return data[(nx * N * N + ny * N + nz) * 27 + i * 9 + j * 3 + k];
    }

    const double& operator()(int nx, int ny, int nz, int i, int j, int k) const {
        return data[(nx * N * N + ny * N + nz) * 27 + i * 9 + j * 3 + k];
    }
};

// Spatial slice class containing all relevant fields
class SpatialSlice {
public:
    ScalarField phi;
    ScalarField K;

    ScalarField alpha;
    VectorField beta;
    VectorField B;

    TensorField gammaTilde;
    TensorField gammaTildeInv;

    TensorField gamma;
    TensorField gammaInv;
    
    TensorField ATilde;
    TensorField ATildeInv;

    Tensor3Field Christoffel;
    Tensor3Field ChristoffelTilde;
    TensorField Ricci;

    VectorField GammaTilde;

    TensorField aux;

    // The main fields are: phi, K, gammaTilde, ATilde, GammaTilde

    SpatialSlice& operator+=(const SpatialSlice& rhs) {
        phi += rhs.phi;
        K += rhs.K;
        gammaTilde += rhs.gammaTilde;
        ATilde += rhs.ATilde;
        GammaTilde += rhs.GammaTilde;
        return *this;
    }

    SpatialSlice& operator*=(const double& rhs) {
        phi *= rhs;
        K *= rhs;
        gammaTilde *= rhs;
        ATilde *= rhs;
        GammaTilde *= rhs;
        return *this;
    }

    SpatialSlice& operator=(const SpatialSlice& rhs) {
        phi = rhs.phi;
        K = rhs.K;
        gammaTilde = rhs.gammaTilde;
        ATilde = rhs.ATilde;
        GammaTilde = rhs.GammaTilde;
        return *this;
    }
};

#endif // FIELD_HPP
