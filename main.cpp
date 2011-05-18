/*
    main.cpp (Multigrid)
    Jess Robertson, 2011-01-28    
    
    Test code for multigrid library.
*/          
 
#include <iostream>       
#include "multigrid.hpp"

using namespace mgrid;     

void check_boundary_conditions(const int nx, const int nz) {
    // Define a boundary point for neumann and dirichlet conditions
    BoundaryPoint neumannCondition = { neumann, 0.0 }; 
    BoundaryPoint dirichletCondition = { dirichlet, 0.0 }; 
    
    // Make a stack with default settings, set some basic conditions
    Settings settings;
    Stack stack(settings);  
    stack.boundaryConditions.set(leftBoundary, neumannCondition);
    stack.boundaryConditions.set(rightBoundary, dirichletCondition);
    stack.boundaryConditions.set(topBoundary, neumannCondition);   
    
    // Make up some condition for the bottom boundary
    Boundary boundary(nx);
    Boundary::iterator bIter;
    int val = 0;
    for (bIter = boundary.begin(); bIter != boundary.end(); bIter++, val++) {
        bIter->value = val;
        if (val > 5) {
            bIter->conditionType = dirichlet;
        } else {
            bIter->conditionType = neumann;
        }
    }  
    stack.boundaryConditions.set(bottomBoundary, boundary);
    
    // Check that boundary conditions have been set
    BoundaryFlag boundaryFlag = bottomBoundary;
    for (Level level=3; level>=0; level--) {  
        std::cout << "Level: " << level << std::endl;
        foreach(BoundaryPoint pt, stack[level].boundaryConditions.get(boundaryFlag)) 
            std::cout << "{" << pt.conditionType << ", " << pt.value << "} ";
        std::cout << std::endl; 
        std::cout.flush();
    }
}

int main (int argc, char const *argv[]) { 
    // Make a stack with default settings, set some basic conditions
    Settings settings;    
    settings.aspectRatio = 4;
    Stack stack(settings);  
    stack.boundaryConditions.set(leftBoundary, zeroNeumannCondition);
    stack.boundaryConditions.set(rightBoundary, zeroDirichletCondition);
    stack.boundaryConditions.set(topBoundary, zeroNeumannCondition); 
    stack.boundaryConditions.set(bottomBoundary, zeroDirichletCondition);
    
    for (int i=stack.finestLevel; i > stack.coarsestLevel; i--) {
        std::cout << "(" << stack[i].length(0) 
                  << ", " << stack[i].length(1) << ")" << std::endl;
    }
    
    return 0;
}