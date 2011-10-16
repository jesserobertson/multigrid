/*
    mosolov_settings.cpp (Multigrid)
    Jess Robertson, 2010-11-16

    Default settings for multigrid solvers
*/

#include "mosolov_settings.hpp"

// Default settings for Settings  
static const double           defaultAspectRatio          = 2;   
static const int              defaultMaximumIterations    = 400;  
static const mgrid::CycleType defaultMgCycleType          = mgrid::wCycle; 
static const int              defaultMinimimumResolution  = 4;   
static const int              defaultNumberOfGrids        = 8;     
static const unsigned long    defaultPreMGRelaxIter       = 1;
static const unsigned long    defaultPostMGRelaxIter      = 2; 
static const double           defaultResidualTolerance    = 1e-10;
                                                          
// Default settings for MosolovSettings                   
static const double           defaultBinghamNumber        = 0.10;   
static const double           defaultCrustStrength        = 0.15;
static const unsigned long    defaultMaxLagrangeIteration = 1000;  
static const double           defaultLagrangeTolerance    = 1e-6; 
static const double           defaultAugmentingParameter  = 1;  

// Default settings on construction
MosolovSettings::MosolovSettings(): 
    binghamNumber(defaultBinghamNumber),     
    crustStrength(defaultCrustStrength),
    maxLagrangeIteration(defaultMaxLagrangeIteration),
    augmentingParameter(defaultAugmentingParameter),
    lagrangeTolerance(defaultLagrangeTolerance) 
{
    multigridSettings.aspectRatio       = defaultAspectRatio; 
    multigridSettings.numberOfGrids     = defaultNumberOfGrids; 
    multigridSettings.minimumResolution = defaultMinimimumResolution;
    multigridSettings.residualTolerance = defaultResidualTolerance;
    multigridSettings.maximumIterations = defaultMaximumIterations;
    multigridSettings.mgCycleType       = defaultMgCycleType;
    multigridSettings.preMGRelaxIter    = defaultPreMGRelaxIter;
    multigridSettings.postMGRelaxIter   = defaultPostMGRelaxIter; 
}