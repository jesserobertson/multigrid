/*
    multigrid_linear.hpp (Multigrid)
    Jess Robertson, 2011-01-28
    
    Linear multigrid solver class
*/                               

#ifndef MULTIGRID_LINEAR_HPP_M1Z9C6EC
#define MULTIGRID_LINEAR_HPP_M1Z9C6EC

#include "multigrid_base.hpp"

namespace mgrid {  

class LinearMultigrid: public MultigridBase {
public:
    LinearMultigrid(const Settings& settings):
        MultigridBase::MultigridBase(settings) {};
    virtual ~LinearMultigrid() {};   

    // Multigrid method
    virtual void multigrid();          
};

} // end namespace mgrid

#endif /* end of include guard: MULTIGRID_LINEAR_HPP_M1Z9C6EC */
