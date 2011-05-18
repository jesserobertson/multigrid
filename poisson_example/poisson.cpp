/*
    poisson.cpp (poisson)
    Jess Robertson, 2010-11-17

    Write methods for Poisson class
*/                                 

#include "poisson.hpp"    

using namespace mgrid;

Poisson::Poisson(const Settings& settings): LinearMultigrid::LinearMultigrid(settings) {
    nxfine = solution[finestLevel].rows(); 
    nzfine = solution[finestLevel].columns();
    std::cout << " -- Initial dimensions: (" << nxfine << ", " << nzfine 
              << "), with aspect: " << settings.aspectRatio << std::endl;
    
    // Set boundary conditions for velocity array 
    solution.boundaryConditions.set(leftBoundary,   zeroNeumannCondition);
    solution.boundaryConditions.set(rightBoundary,  zeroDirichletCondition);
    solution.boundaryConditions.set(topBoundary,    zeroNeumannCondition);    
    solution.boundaryConditions.set(bottomBoundary, zeroDirichletCondition); 
    
    // Set source  
    source[finestLevel] = -1.0; sourceIsSet = true;    
}    
Poisson::~Poisson() { /* pass */ }

std::string Poisson::filename(std::string root) {
    std::ostringstream name;  
    name.precision(1);  // Print variables to one decimal place
    name << root << "A" << std::fixed << aspect; 
    return name.str();
} 

void Poisson::solve() {       
    // Solve using linear multigrid method
    multigrid(); 
}