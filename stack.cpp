/*
    stack.cpp (Multigrid)
    Jess Robertson, 2010-09-15   
    
    Implementation of Stack class
*/          

#include "stack.hpp"                  
                    
// Ctor
mgrid::Stack::Stack(const mgrid::Settings& s): 
    finestLevel(s.numberOfGrids-1), aspect(s.aspectRatio), 
    nGrids(s.numberOfGrids), minRes(s.minimumResolution)
{
    // Generate grid stack
    resize(nGrids); 
                          
    // Generate coarsest level
    int nx, nz;
    if (aspect <= 1.99999999999999999999999) {
        nx = minRes;
        nz = int(round(2*(minRes-1)/aspect) + 1);
    } else {
        nx = int(round(aspect*(minRes-1)/2.0 + 1));
        nz = minRes;
    }   
    (*this)[coarsestLevel].resize(aspect, nx, nz); 
    (*this)[coarsestLevel] = 0;

    // Proceed for all other grid sizes, caching geometry as we go...
    for (Level level=1; level<=finestLevel; level++) {   
        nx = 2*(nx - 1) + 1; 
        nz = 2*(nz - 1) + 1;
        (*this)[level].resize(aspect, nx, nz);      
        (*this)[level] = 0;
    }     
    
    // Set default boundary conditions
    boundaryConditions.apply_default_conditions(nx, nz);
    
    // Make boundary conditions references of boundary conditions on finest level  
    foreach(BoundaryFlag boundaryFlag, allBoundaryFlags)
        _update_boundary_conditions(boundaryFlag);
}

void mgrid::Stack::_update_boundary_conditions(BoundaryFlag boundaryFlag) {
    // Updates boundary conditions on given boundary at all levels
    // Get length of finest boundary condition    
    int fineLength = boundaryConditions.get(boundaryFlag).extent(0);   
    
    // Top level is a direct reference
    (*this)[finestLevel].boundaryConditions.get(boundaryFlag)\
        .reference(boundaryConditions.get(boundaryFlag));
    
    // Make references for all other levels. This means that this method should
    // probably only be called once, on construction of the Stack instance, as changing
    // the boundary conditions should be reflected in all the other levels 
    int strideLength = 1;
    for (Level level=finestLevel-1; level>=coarsestLevel; level--) {
        strideLength *= 2; 
        blitz::Range strideDomain(0,fineLength,strideLength);   
        (*this)[level].boundaryConditions.get(boundaryFlag)\
            .reference(boundaryConditions.get(boundaryFlag)(strideDomain));
    }   
}