/*
    multigrid_nonlinear.hpp (Multigrid)
    Jess Robertson, 2011-01-28
    
    Linear multigrid solver class
*/                               

#ifndef MULTIGRID_NONLINEAR_HPP_J1PG82P8
#define MULTIGRID_NONLINEAR_HPP_J1PG82P8

#include "multigrid_base.hpp"

namespace mgrid {  

class NonlinearMultigrid: public MultigridBase {
public:
    NonlinearMultigrid(const Settings& settings):
        MultigridBase::MultigridBase(settings),
        truncError(settings),
        rightHandSide(settings) {};
    virtual ~NonlinearMultigrid () {};    
    
    // Multigrid method
    virtual void multigrid();        
    
protected:
    Stack truncError;    // Truncation error    
    Stack rightHandSide; // Right hand side of different equations
}; 

} // end namespace mgrid

#endif /* end of include guard: MULTIGRID_NONLINEAR_HPP_J1PG82P8 */