#ifndef FIELDDATA_HPP
#define FIELDDATA_HPP

#include "fieldTypes.hpp"

class FieldData {
public:
    // Get the singleton instance
    static FieldData& getInstance();

    // Accessors for the spatial slices
    SpatialSlice& getCurrentSlice();
    SpatialSlice& getFutureSlice();
    SpatialSlice& getAuxSlice();

    // Delete copy constructor and assignment operator
    FieldData(const FieldData&) = delete;
    FieldData& operator=(const FieldData&) = delete;

private:
    // Private constructor
    FieldData();

    // The three spatial slices
    SpatialSlice current;
    SpatialSlice future;
    SpatialSlice aux;
};

#endif // FIELDDATA_HPP
