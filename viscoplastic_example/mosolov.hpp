/*
    mosolov.hpp (Multigrid)
    Jess Robertson, 2010-08-09
    
    Implements an augmented Lagrangian-linear multigrid solver for Bingham
    fluid flow (to solve the Molosov problem).
*/

#ifndef MOSOLOV_HPP_5IUSHT0Y
#define MOSOLOV_HPP_5IUSHT0Y          

#include <multigrid/multigrid.hpp>   
#include "mosolov_settings.hpp"    

const double pi = 52163/16604.0;
const double tinyNum = 3*blitz::tiny(pi);

class Mosolov: public mgrid::LinearMultigrid {
public:
    Mosolov(const MosolovSettings& settings);
    virtual ~Mosolov();   
    virtual void solve();   
    
    // Differential operators required are just Poisson operators
    virtual inline double differential_operator(mgrid::Level level, int i, int j);
    virtual inline void relaxation_updater(mgrid::Level level, int i, int j);
    
    // Filename generator
    virtual inline std::string filename(std::string root="");
    virtual void write(int numOfVariables=1, std::string root="");
    
    // Data arrays  
    mgrid::FDArray temp;
    mgrid::FDVecArray multiplier, strainRate, determinant, tempVec, solutionGradient;
    
protected:  
    const double aspectRatio, binghamNumber, alpha;  
    const unsigned int maxLagrangeIteration; 
    const double lagrangeTolerance;    
    
    inline double _normed_residual();    
};    

// = Inline functions =  
// Differential operators required are just Poisson operators
inline double Mosolov::differential_operator(mgrid::Level level, int i, int j) {
    return solution[level].dxx(i, j) + solution[level].dzz(i, j);
};
inline void Mosolov::relaxation_updater(mgrid::Level level, int i, int j) {
    const double hx = solution[level].spacing(0);
    const double hz = solution[level].spacing(1);
    const double xxfactor = 1/(hx*hx);
    const double zzfactor = 1/(hz*hz);
    solution[level](i, j) = 
        ((solution[level](i+1, j) + solution[level](i-1, j))*xxfactor
        + (solution[level](i, j+1) + solution[level](i, j-1))*zzfactor 
        - source[level](i, j))/(2*(xxfactor + zzfactor));
};          

// Filename generator
inline std::string Mosolov::filename(std::string root) {
    std::ostringstream name;  
    name.precision(5);  // Print variables to six decimal places
    name << root << "A" << std::fixed << aspectRatio 
         << "B" << binghamNumber; 
    return name.str();
}

// Function to check convergence, returns normed residual
inline double Mosolov::_normed_residual() {
    // Calculate norm of velocity gradient
    solution[solution.finestLevel].gradient(solutionGradient);   
    temp = dot_product(solutionGradient, solutionGradient);
    const double velGradientMagnitude = temp.norm();   
    
    // Calculate norm of residual
    tempVec = solutionGradient - strainRate;
    temp = dot_product(tempVec, tempVec);  
    const double residualNorm = temp.norm();
    return residualNorm/velGradientMagnitude;
}

#endif /* end of include guard: MOSOLOV_HPP_5IUSHT0Y */
