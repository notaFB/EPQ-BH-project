
#ifndef INITIAL_CONDITIONS_HPP
#define INITIAL_CONDITIONS_HPP

#include "fieldData.hpp"

class InitialConditions {
public:
    static void setgammaTilde(SpatialSlice& slice);
    static void setATilde(SpatialSlice& slice);
    static void setK(SpatialSlice& slice);
    static void setGammaTilde(SpatialSlice& slice);
    static void setPhi(SpatialSlice& slice);
    static void setalpha(SpatialSlice& slice);
    static void setBeta(SpatialSlice& slice);
    static void setB(SpatialSlice& slice);

    static void setInitialConditions(SpatialSlice& slice);

};



#endif