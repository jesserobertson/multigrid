/*
    settings.cpp (Multigrid)
    Jess Robertson, 2011-01-28
*/          

#include "settings.hpp"

// Default settings for Settings  
static const double           defaultAspectRatio             = 1;   
static const int              defaultMaximumIterations       = 400;  
static const mgrid::CycleType defaultMgCycleType             = mgrid::wCycle; 
static const int              defaultMinimimumResolution     = 4;   
static const int              defaultNumberOfGrids           = 8;     
static const unsigned long    defaultPreMGRelaxIter          = 1;
static const unsigned long    defaultPostMGRelaxIter         = 2; 
static const double           defaultResidualTolerance       = 1e-10;

// Apply default settings on construction
mgrid::Settings::Settings():
    aspectRatio(defaultAspectRatio), 
    numberOfGrids(defaultNumberOfGrids), 
    minimumResolution(defaultMinimimumResolution),
    residualTolerance(defaultResidualTolerance),
    maximumIterations(defaultMaximumIterations),
    mgCycleType(defaultMgCycleType),
    preMGRelaxIter(defaultPreMGRelaxIter),
    postMGRelaxIter(defaultPostMGRelaxIter) { /* pass */ }