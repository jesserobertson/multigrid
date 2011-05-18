/*
    multigrid_base.hpp (Multigrid)
    Jess Robertson, 2011-01-28
*/             

#ifndef MULTIGRID_BASE_HPP_TVC215N7
#define MULTIGRID_BASE_HPP_TVC215N7

#include <netcdfcpp.h> 

#include "types.hpp" 
#include "multigrid_exceptions.hpp"
#include "utilities.hpp"
#include "fdarray.hpp"
#include "fdvecarray.hpp"
#include "stack.hpp" 
#include "settings.hpp"

namespace mgrid {

class MultigridBase {
public:         
    MultigridBase(const Settings& settings);       
    virtual ~MultigridBase() {}    
    
    // Setters and getters 
    inline FDArray& get_result(); 
    template <typename T> inline void initial_guess(T arg); 
    inline FDArray& source_term(); 
    template <typename T> inline void source_term(T arg);
    
    // Evaluation methods 
    inline void evaluate_operator(Level level, FDArray& result); 
    inline void evaluate_residual(Level level, FDArray& result);                     
    
    // Relaxation methods  
    void relax(const Level level, const unsigned long N);
    void relax(const Level level, const double tolerance);        
    
    // Multigrid solver method, overwritten by LinearMultigrid and 
    // NonlinearMultigrid classes, and solve method which should be 
    // overwritten by subclasses of Linear- and NonlinearMultigrid if
    // different behavior than just calling multigrid() is desired. 
    virtual inline void multigrid() { /* pass */ }  
    virtual inline void solve() { multigrid(); }        
    
    // Other overloaded methods
    virtual double differential_operator(Level, int, int)=0;
    virtual void relaxation_updater(Level, int, int)=0;
    
    // Writing methods & file name generator
    virtual void write(int numOfVariables, std::string root="");   
    virtual std::string filename(std::string root="")=0;
  
protected:      
    // Data  
    Stack solution, source;         // Grids for solution and source term
    Stack temp;                     // Extra storage for multigrid solver
    const int cycleType;            // Type of FMG-cycling used
    const unsigned long preRelax;   // Num of pre-corection relaxations to use
    const unsigned long postRelax;  // Num of post-corection relaxations to use    
    const double residualTolerance; // For convergence testing
    const double maxIterations;     // Maxium number of iterations allowed    
    const double aspect;            // Aspect ratio  
    int finestLevel, coarsestLevel, nxfine, nzfine;  // Grid geometry        
    bool sourceIsSet;               // Has the source term been provided?
    bool initialIsSet;              // Has an initial value for the solution
                                    // been provided? (This can be useful for 
                                    // solve routines, which may generate 
                                    // their own initial values otherwise).
    
private: 
    double residualSum, normSum;  
    Deriv du;
};    

// Setters and getters           
inline FDArray& MultigridBase::get_result() {
    return solution[finestLevel];
}
template <typename T> inline void MultigridBase::initial_guess(T arg) {
    solution[finestLevel] = arg;
    initialIsSet = true;
} 
inline FDArray& MultigridBase::source_term() {
    return source[finestLevel];
} 
template <typename T> inline void MultigridBase::source_term(T arg) {
    source[finestLevel] = arg;  
    sourceIsSet = true;
}

// Evaluation methods
inline void MultigridBase::evaluate_operator(Level level, FDArray& result) {
    ARRAY_LOOP(result) 
        result(i, j) = differential_operator(level, i, j); 
}   
inline void MultigridBase::evaluate_residual(Level level, FDArray& result) {
    ARRAY_LOOP(result)  
        result(i, j) = source[level](i, j) 
            - differential_operator(level, i, j); 
}   

} // end namespace mgrid

#endif /* end of include guard: MULTIGRID_BASE_HPP_TVC215N7 */
