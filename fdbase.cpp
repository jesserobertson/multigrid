/*
    fdbase.cpp (Multigrid)
    Jess Robertson, 2011-01-28
*/

#include "fdbase.hpp"

mgrid::FDBase::FDBase() { /* pass */ }

void mgrid::FDBase::calculate_geometry(const double aspectRatio, const int nx, const int nz) {
    // Copy over data
    (*this).aspectRatio = aspectRatio; 
    (*this).nx = nx; 
    (*this).nz = nz;                      
    
    // Calculate spacings
    hx = aspectRatio/double(nx-1);      
    hz = 1/double(nz-1); 
    
    // Pre-calculate factors for derivatives   
    xfactor = 1.0/(2*hx);       zfactor = 1.0/(2*hz);
    xxfactor = 1.0/(hx*hx);     zzfactor = 1.0/(hz*hz);
    xzfactor = 1.0/(4*hx*hz);         
}   
