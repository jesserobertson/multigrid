/*
    boundary_conditions.hpp (Multigrid)
    Jess Robertson, 2010-09-16
    
    Interface for boundary conditions classes
*/                       

#ifndef BOUNDARY_CONDITIONS_HPP_3I5KIZNP
#define BOUNDARY_CONDITIONS_HPP_3I5KIZNP

#include "types.hpp" 
#include "multigrid_exceptions.hpp"
#include "utilities.hpp"   

namespace mgrid {

// = BoundaryConditions class interface =
class BoundaryConditions {
public:
    // Ctors
    BoundaryConditions();
    virtual ~BoundaryConditions();
    BoundaryConditions(const BoundaryConditions& copyFrom); 
    const BoundaryConditions& operator=(const BoundaryConditions& copyFrom);
    
    // Settors                       
    inline Boundary& get(BoundaryFlag boundaryFlag) {
        return boundaries(boundaryFlag);
    }                                           
    void set(BoundaryFlag boundaryFlag, const Boundary& referent);
    void set(BoundaryFlag boundaryFlag, const BoundaryPoint& pt);
    inline void resize(const int nx, const int nz) {
        boundaries(topBoundary).resize(nx);
        boundaries(bottomBoundary).resize(nx);
        boundaries(leftBoundary).resize(nz);
        boundaries(rightBoundary).resize(nz);
    } 
    
    // Apply default boundary conditions  
    void apply_default_conditions(const int nx, const int nz);       
    
    // Reference another set of BoundaryConditions
    void reference(const BoundaryConditions& referent);
    
private:
    blitz::TinyVector<Boundary, 4> boundaries;
};

// = Some point definitions whch get used a lot
const BoundaryPoint zeroNeumannCondition = { neumann, 0.0 }; 
const BoundaryPoint zeroDirichletCondition = { dirichlet, 0.0 };

} // end namespace Multigrid  

#endif /* end of include guard: BOUNDARY_CONDITIONS_HPP_3I5KIZNP */

