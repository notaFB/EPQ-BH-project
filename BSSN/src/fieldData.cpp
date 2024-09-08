#include "fieldData.hpp"

// Initialize the singleton instance
FieldData& FieldData::getInstance() {
    static FieldData instance;
    return instance;
}

// Constructor for FieldData, initializes slices
FieldData::FieldData()
    : v1(), k1(), k2(), k3(), k4(), aux() {}

// Accessors for each slice
SpatialSlice& FieldData::getV1Slice() {
    return v1;
}

SpatialSlice& FieldData::getk1Slice() {
    return k1;
}

SpatialSlice& FieldData::getk2Slice() {
    return k2;
}

SpatialSlice& FieldData::getk3Slice() {
    return k3;
}

SpatialSlice& FieldData::getk4Slice() {
    return k4;
}

SpatialSlice& FieldData::getAuxSlice() {
    return aux;
}

