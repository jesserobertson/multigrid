/*
    poisson.hpp (poisson)
    Jess Robertson, 2010-08-08
*/

#ifndef POISSON_HPP_4E5A9H5B
#define POISSON_HPP_4E5A9H5B    

#include <multigrid/multigrid.hpp>

// = Poisson class interface =
class Poisson: public mgrid::LinearMultigrid {
public:
    Poisson(const mgrid::Settings& settings);  
    virtual ~Poisson(); 
    
    // Solution routines    
    virtual void solve();
    virtual inline double differential_operator(mgrid::Level level, int i, int j); 
    virtual inline void relaxation_updater(mgrid::Level level, int i, int j);   
    
    // Filename generator
    virtual std::string filename(std::string root="");  
}; 

// = Inline functions =     
// Differential operators
inline double Poisson::differential_operator(mgrid::Level level, int i, int j) {
    return solution[level].dxx(i, j) + solution[level].dzz(i, j);
};
inline void Poisson::relaxation_updater(mgrid::Level level, int i, int j) { 
    const double hx = solution[level].spacing(0);
    const double hz = solution[level].spacing(1);
    const double xxfactor = 1/(hx*hx);
    const double zzfactor = 1/(hz*hz);
    solution[level](i, j) = 
        ((solution[level](i+1, j) + solution[level](i-1, j))*xxfactor
        + (solution[level](i, j+1) + solution[level](i, j-1))*zzfactor 
        - source[level](i, j))/(2*(xxfactor + zzfactor));
};  

#endif /* end of include guard: POISSON_HPP_4E5A9H5B */