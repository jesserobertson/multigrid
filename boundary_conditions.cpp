/*
    boundary_conditions.cpp (Multigrid)
    Jess Robertson, 2011-01-28
    
    Implementation of a simple boundary conditions class
*/ 

#include "boundary_conditions.hpp"    

// Ctors etc
mgrid::BoundaryConditions::BoundaryConditions() { /* pass */ }
mgrid::BoundaryConditions::~BoundaryConditions() { /* pass */ }
mgrid::BoundaryConditions::BoundaryConditions(const BoundaryConditions& copyFrom) { /* pass */ } 
const mgrid::BoundaryConditions& 
    mgrid::BoundaryConditions::operator=(const BoundaryConditions& copyFrom) { /* pass */ }

// Accessor methods 
void mgrid::BoundaryConditions::set(BoundaryFlag bFlag, const Boundary& referent) {
    // Check lengths before assignment
    if (boundaries(bFlag).extent(0) != referent.extent(0)) {
        throw InvalidBoundaryCondition();
    } else { 
        Boundary::iterator boundIter = boundaries(bFlag).begin();  
        for (Boundary::const_iterator refIter = referent.begin(); 
             refIter != referent.end(); 
             refIter++, boundIter++) 
        {
            boundIter->value = refIter->value;
            boundIter->conditionType = refIter->conditionType;
        }       
    }
} 
void mgrid::BoundaryConditions::set(BoundaryFlag bFlag, const BoundaryPoint& pt) {
    boundaries(bFlag) = pt;      
}      

// Some default boundary conditions
void mgrid::BoundaryConditions::apply_default_conditions(const int nx, const int nz) {
    // Set some simple boundary points
    BoundaryPoint neumannConditionPt = { neumann, 0.0 }; 
    BoundaryPoint dirichletConditionPt = { dirichlet, 0.0 };
    
    // Initialise with default boundary conditions   
    resize(nx, nz);
    set(leftBoundary, neumannConditionPt);
    set(rightBoundary, dirichletConditionPt);
    set(topBoundary, neumannConditionPt);
    set(bottomBoundary, dirichletConditionPt);
}

// Reference another set of BoundaryConditions
void mgrid::BoundaryConditions::reference(const BoundaryConditions& referent) {
    foreach(BoundaryFlag boundaryFlag, allBoundaryFlags)
        boundaries(boundaryFlag).reference(referent.boundaries(boundaryFlag));
}