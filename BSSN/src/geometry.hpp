#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "fieldTypes.hpp"
#include "utilities.hpp"


class Geometry {
public:

static void computeInverseMetric(TensorField& gamma, TensorField& gammaInv);
static void computeInverseTensor(TensorField& T, TensorField& TInv, TensorField& gamma);
//Christoffel symbols of the second kind
static void computeChristoffel(TensorField& gamma, TensorField& gammaInv, Tensor3Field& Christoffel);

static void computerRicci(SpatialSlice& slice);
static void computegamma(SpatialSlice& slice);

static void computeGeometry(SpatialSlice& slice);

};


#endif