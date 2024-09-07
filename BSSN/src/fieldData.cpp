#include "fieldData.hpp"

// Initialize the singleton instance
FieldData& FieldData::getInstance() {
    static FieldData instance;
    return instance;
}

// Constructor for FieldData, initializes slices
FieldData::FieldData()
    : v1(), v2(), v3(), v4() {}

// Accessors for each slice
SpatialSlice& FieldData::getV1Slice() {
    return v1;
}

SpatialSlice& FieldData::getV2Slice() {
    return v2;
}

SpatialSlice& FieldData::getV3Slice() {
    return v3;
}

SpatialSlice& FieldData::getV4Slice() {
    return v4;
}

