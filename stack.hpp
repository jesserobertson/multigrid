/*
    stack.hpp (Multigrid)
    Jess Robertson, 2010-09-15     
    
    Simple grid stack class for multigrid problems.
*/                            

#ifndef STACK_HPP_1FKL9HHJ
#define STACK_HPP_1FKL9HHJ

#include <vector> 
#include <boost/tuple/tuple.hpp>  

#include "types.hpp" 
#include "multigrid_exceptions.hpp"
#include "utilities.hpp"                   
#include "fdarray.hpp"   
#include "settings.hpp"
#include "boundary_conditions.hpp"  

namespace mgrid {

// = Transfer operators =
/*  Interpolation and restriction operators which operate on FDArray instances:
    -- restrictor applies a restriction operator which transfers the 
       data on the current level (as specified by the private variable 
       level) to the next coarsest level using fully weighted 
       restriction. It updates the data on the next coarsest grid in the 
       stack and the value of currentLevel.
    -- interpolator applies an interpolation operator which transfers 
       the data on the current level in the grid stack to the next finest 
       level using bilinear interpolation. It updates the data on the next 
       finest grid in the stack and the value of currentLevel.
*/
inline void restriction_operator(FDArray coarse, FDArray fine) {  
    int nxc, nzc, nxf, nzf;   
    boost::tie(nxc, nzc) = to_tuple(coarse.shape());
    boost::tie(nxf, nzf) = to_tuple(fine.shape());   
    
    // Perform restriction over center of grid  
    blitz::Range fi(2, nxf-3, 2); 
    blitz::Range fj(2, nzf-3, 2); 
    blitz::Range ci(1, nxc-2); 
    blitz::Range cj(1, nzc-2);
    coarse(ci, cj) = (4*(fine(fi, fj))
        + 2*(fine(fi+1,fj) + fine(fi-1,fj)+ fine(fi,fj+1) + fine(fi,fj-1))
        + 1*(fine(fi+1,fj+1) + fine(fi+1,fj-1) + fine(fi-1,fj+1) 
            + fine(fi-1,fj-1)))/16.0; 
    
    // Perform restriction at boundaries
    coarse(ci, 0) = (4*fine(fi, 0)
        + 2*(fine(fi-1, 0) + fine(fi+1, 0) + fine(fi, 1))
        + 1*(fine(fi-1, 1) + fine(fi+1, 1)))/12.0;
    coarse(ci, nzc-1) = (4*fine(fi, nzf-1)
        + 2*(fine(fi-1, nzf-1) + fine(fi+1, nzf-1) + fine(fi, nzf-2)) 
        + 1*(fine(fi-1, nzf-2) + fine(fi+1, nzf-2)))/12.0;
    coarse(0, cj) = (4*fine(0, fj)
        + 2*(fine(0, fj-1) + fine(0, fj+1) + fine(1, fj))
        + 1*(fine(1, fj-1) + fine(1, fj+1)))/12.0;
    coarse(nxc-1, cj) = (4*fine(nxf-1, fj)
        + 2*(fine(nxf-1, fj-1) + fine(nxf-1, fj+1) + fine(nxf-2, fj))
        + 1*(fine(nxf-2, fj-1) + fine(nxf-2, fj+1)))/12.0; 
    
    // Perform restriction at corners
    coarse(0,0) = (4*fine(0,0) + 2*(fine(1,0) + fine(0,1)) + 1*fine(1,1))/9.0;
    coarse(nxc-1,0) = (4*fine(nxf-1,0) + 2*(fine(nxf-2,0) + fine(nxf-1,1))
        + 1*fine(nxf-2,1))/9.0;
    coarse(0,nzc-1) = (4*fine(0,nzf-1) + 2*(fine(1,nzf-1) + fine(0,nzf-2)) 
        + 1*fine(1,nzf-2))/9.0;
    coarse(nxc-1,nzc-1) = (4*fine(nxf-1,nzf-1) + 2*(fine(nxf-2,nzf-1) 
        + fine(nxf-1,nzf-2)) + 1*fine(nxf-2,nzf-2))/9.0;
}
inline void interpolation_operator(FDArray coarse, FDArray fine) {
    int nxc, nzc, nxf, nzf;   
    boost::tie(nxc, nzc) = to_tuple(coarse.shape());
    boost::tie(nxf, nzf) = to_tuple(fine.shape());  
    blitz::Range i(0, nxf-1, 2);
    blitz::Range j(0, nzf-1, 2); 
    blitz::Range m(1, nxf-2, 2);
    blitz::Range n(1, nzf-2, 2);
            
    // Copy over data directly 
    for (int ii=0; ii<nxc; ii++) 
        for (int jj=0; jj<nzc; jj++)
            fine(2*ii, 2*jj) = coarse(ii, jj);
    
    // Interpolation over center of grid
    fine(i, n) = 0.5*(fine(i, n-1) + fine(i, n+1));
    fine(m, j) = 0.5*(fine(m-1, j) + fine(m+1, j));
    fine(m, n) = 0.25*(fine(m+1, n+1) + fine(m+1, n-1) + fine(m-1, n+1) 
        + fine(m-1, n-1));     
        
    // Interpolation at boundaries of grid  
    fine(0,n) = 0.5*(fine(0, n-1) + fine(0, n+1));
    fine(nxf-1,n) = 0.5*(fine(nxf-1, n-1) + fine(nxf-1, n+1));
    fine(m, 0) = 0.5*(fine(m-1, 0) + fine(m+1, 0));
    fine(m, nzf-1) = 0.5*(fine(m-1, nzf-1) + fine(m+1, nzf-1));   
}
    
// = Stack class interface =
class Stack: public std::vector<mgrid::FDArray> {
public:
    Stack(const Settings& s);  
    virtual ~Stack () {};          

    // Some useful attributes
    const Level finestLevel;              // level of finest grid  
    static const Level coarsestLevel = 0; // level of coarsest grid      

    // Methods
    void write(std::string fileString); 
    inline void coarsen(Level level);
    inline void coarsen(Level level, FDArray& result);
    inline void refine(Level level);
    inline void refine(Level level, FDArray& result); 
    
    // Boundary conditions methods      
    BoundaryConditions boundaryConditions;
    void _update_boundary_conditions(BoundaryFlag boundaryFlag);  

private:
    // Boundary conditions data  
    // void _update_boundary_conditions(BoundaryFlag boundaryFlag); 
    
    // Geometry data              
    const double aspect;          // aspect ratio of grid domain  
    const int nGrids;             // number of grid levels required
    int minRes;                   // minimum resolution parameter        
};                   

// = Inline methods for Stack class =     
inline void Stack::coarsen(Level level) {
    restriction_operator((*this)[level - 1], (*this)[level]);
}
inline void Stack::coarsen(Level level, FDArray& result) {
    restriction_operator(result, (*this)[level]);
}
inline void Stack::refine(Level level) {
    interpolation_operator((*this)[level], (*this)[level + 1]);
}        
inline void Stack::refine(Level level, FDArray& result) {
    interpolation_operator((*this)[level], result);
}
       
} // end namespace multigrid        

#endif /* end of include guard: STACK_HPP_1FKL9HHJ */
