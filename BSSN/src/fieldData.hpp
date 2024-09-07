#ifndef FIELDDATA_HPP
#define FIELDDATA_HPP

#include "fieldTypes.hpp"

class FieldData {
public:
    // Get the singleton instance
    static FieldData& getInstance();

    // Accessors for the four spatial slices
    SpatialSlice& getV1Slice();
    SpatialSlice& getV2Slice();
    SpatialSlice& getV3Slice();
    SpatialSlice& getV4Slice();

    // Delete copy constructor and assignment operator
    FieldData(const FieldData&) = delete;
    FieldData& operator=(const FieldData&) = delete;

private:
    // Private constructor
    FieldData();

    // The four spatial slices
    SpatialSlice v1;
    SpatialSlice v2;
    SpatialSlice v3;
    SpatialSlice v4;
};

#endif // FIELDDATA_HPP
