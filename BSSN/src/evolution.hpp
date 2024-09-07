#ifndef EVOLUTION_HPP
#define EVOLUTION_HPP


#include"fieldData.hpp"

// dS/dt = F(S)

class Evolution {
public:
    Evolution(){}

    static void F(SpatialSlice& slice, SpatialSlice& sliceOut) {
        
        F_phi(slice, sliceOut);
        F_K(slice, sliceOut);
        F_alpha(slice, sliceOut);
        F_beta(slice, sliceOut);
        F_B(slice, sliceOut);
        F_gammaTilde(slice, sliceOut);
        F_ATilde(slice, sliceOut);
        F_GammaTilde(slice, sliceOut);

    }


    static void F_phi(SpatialSlice& slice, SpatialSlice& sliceOut);
    static void F_K(SpatialSlice& slice, SpatialSlice& sliceOut);
    static void F_alpha(SpatialSlice& slice, SpatialSlice& sliceOut);
    static void F_beta(SpatialSlice& slice, SpatialSlice& sliceOut);
    static void F_B(SpatialSlice& slice, SpatialSlice& sliceOut);
    static void F_gammaTilde(SpatialSlice& slice, SpatialSlice& sliceOut);
    static void F_ATilde(SpatialSlice& slice, SpatialSlice& sliceOut);
    static void F_GammaTilde(SpatialSlice& slice, SpatialSlice& sliceOut);

    static void RK4_step();

};



#endif