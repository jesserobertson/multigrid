/*
    mosolov_settings.hpp (Multigrid)
    Jess Robertson, 2010-09-19
    
    Argument struct to save passing around settings information
*/                                                      

#ifndef MOSOLOV_SETTINGS_HPP_I32RWBU4
#define MOSOLOV_SETTINGS_HPP_I32RWBU4

#include <multigrid/multigrid.hpp>

// Settings struct
struct MosolovSettings {
    // Linear Multigrid settings
    mgrid::Settings multigridSettings;
        
    // Settings specific to Mosolov class
    double binghamNumber;    
    double crustStrength;
    double augmentingParameter;
    unsigned long maxLagrangeIteration;   
    double lagrangeTolerance;
        
    // Set default values (in .cpp file) on construction 
    MosolovSettings();
};    

#endif /* end of include guard: MOSOLOV_SETTINGS_HPP_I32RWBU4 */
