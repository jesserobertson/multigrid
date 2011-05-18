/*
    settings.hpp (Multigrid)
    Jess Robertson, 2010-09-19
    
    Argument struct to save passing around settings information
*/                                                      

#ifndef SETTINGS_HPP_7RZE6Z3D
#define SETTINGS_HPP_7RZE6Z3D       

#include "types.hpp" 
#include "multigrid_exceptions.hpp"

namespace mgrid {
    
// Settings struct
struct Settings { 
    double aspectRatio; 
    int numberOfGrids; 
    int minimumResolution;
    double residualTolerance;
    int maximumIterations;
    CycleType mgCycleType;
    unsigned long preMGRelaxIter;
    unsigned long postMGRelaxIter;
        
    // Ctor etc
    Settings(); // Default settings in settings.cpp
};      

} // end namespace mgrid
                         
#endif /* end of include guard: SETTINGS_HPP_7RZE6Z3D */
