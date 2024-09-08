#ifndef FIELDDATA_HPP
#define FIELDDATA_HPP

#include "fieldTypes.hpp"

class FieldData {
public:
    // Get the singleton instance
    static FieldData& getInstance();

    // Accessors for the four spatial slices
    SpatialSlice& getV1Slice();
    SpatialSlice& getk1Slice();
    SpatialSlice& getk2Slice();
    SpatialSlice& getk3Slice();
    SpatialSlice& getk4Slice();
    SpatialSlice& getAuxSlice();

    // Delete copy constructor and assignment operator
    FieldData(const FieldData&) = delete;
    FieldData& operator=(const FieldData&) = delete;

private:
    // Private constructor
    FieldData();

    // The four spatial slices
    SpatialSlice v1;
    SpatialSlice k1;
    SpatialSlice k2;
    SpatialSlice k3;
    SpatialSlice k4;
    SpatialSlice aux;
};

#endif // FIELDDATA_HPP
