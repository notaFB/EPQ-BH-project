#include "fieldData.hpp"

// Initialize the singleton instance
FieldData& FieldData::getInstance() {
    static FieldData instance;
    return instance;
}

// Constructor for FieldData, initializes slices
FieldData::FieldData()
    : current(), future(), aux() {}

// Accessors for each slice
SpatialSlice& FieldData::getCurrentSlice() {
    return current;
}

SpatialSlice& FieldData::getFutureSlice() {
    return future;
}

SpatialSlice& FieldData::getAuxSlice() {
    return aux;
}
