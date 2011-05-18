/*
    types.hpp (Multigrid) 
    Jess Robertson, 2010-07-25
    
    Custom type declarations for multigrid solver.
*/     

#ifndef TYPES_HPP_YJLHHNAL
#define TYPES_HPP_YJLHHNAL

#include <blitz/array.h>
#include <boost/tuple/tuple.hpp>      

namespace mgrid {        

// Iteration types
typedef long unsigned int iter_t;

// Flags and boundary condition specifications   
enum CycleType {vCycle = 1, wCycle = 2, threeCycle = 3};

// Deriv structs
typedef struct {
    double dx, dxu, dz, dzu, dxx, dxxu, dzz, dzzu, dxz, dxzu;   
} Deriv;

// Multigrid paramters
typedef int Level;    

// Boundary condition types
enum ConditionType {dirichlet, neumann}; 
struct BoundaryPoint { 
    ConditionType conditionType; 
    double value; 
};
typedef blitz::Array<BoundaryPoint, 1> Boundary;
    
// Boundary flags
enum BoundaryFlag {leftBoundary, rightBoundary, topBoundary, bottomBoundary};
const blitz::TinyVector<BoundaryFlag, 4> 
    allBoundaryFlags(leftBoundary, rightBoundary, topBoundary, bottomBoundary);   

// Some typedefs to ease polymorphism with array methods
typedef blitz::Array<double, 3>::T_array VecArrayType; 
typedef blitz::Array<double, 2>::T_array ArrayType;

} // end namespace mgrid
    
#endif /* end of include guard: TYPES_HPP_YJLHHNAL */