/*
    multigrid_nonlinear.cpp (Multigrid)
    Jess Robertson, 2011-01-28
*/                            

#include "multigrid_nonlinear.hpp"

void mgrid::NonlinearMultigrid::multigrid() {
    // Check that the source array has been set                          
    if (not(sourceIsSet)) return;
    
    // Constants                          
    const Level finestLevel = solution.finestLevel;
    const Level coarsestLevel = solution.coarsestLevel;
    
    // Initialise initial guess
    for (Level level=finestLevel; level>0; level--) {
        solution.coarsen(level); 
        source.coarsen(level);                         
    }

    // Solve on coarsest level 
    rightHandSide[coarsestLevel] = source[coarsestLevel];    
    relax(coarsestLevel, residualTolerance);        

    // Full Multigrid loop
    for (Level fineLevel=1; fineLevel<=finestLevel; fineLevel++) {
        // V-cycle loop at each (successively finer) level
        solution.refine(fineLevel-1); // interpolate solution to next level   
        rightHandSide[fineLevel] = source[fineLevel];  // set up rRHS
        for (int cycle=0; cycle < cycleType; cycle++) {         
            // Downstroke of cycle:
            //  -- New solution: u(2h) = R.u(h)     
            //  -- Truncation Error: t = L(R.u(h)) - R(L.u(h)) 
            //  -- New rhs: f(2h) = R.f(h) + t
            for (Level level=fineLevel; level>coarsestLevel; level--) {
                // Do pre-correction relaxation on current level 
                relax(level, preRelax);
            
                // Calculate approximate truncation error for current
                // discretisation on the finest grid: t = L(R.u(h)) - R(L.u(h))
                evaluate_operator(level, temp[level]);     // L.u(h)      
                temp.coarsen(level);                       // temp <- R(L.u(h)) 
                solution.coarsen(level);                   // u(2h) <- R.u(h)
                evaluate_operator(level-1, truncError[level-1]); // L(R.u(h))      
                truncError[level-1] -= temp[level-1];      // t
                
                // // Estimate truncation error based on Taylor expansion of
                // // PDE residual (see NRC, pp. ?): |t(h)| ≈ |t(2h)|/3  
                // if (level == fineLevel)    
                //     double truncErrorNorm = (temp[level-1].norm())/3.0;    
                
                // Calculate source: f(2h) = R.f(h) + L(R.u(h)) - R(L.u(h))
                source.coarsen(level);
                source[level-1] += truncError[level-1];   
            }   
    
            // Solve on coarsest level  
            relax(coarsestLevel, residualTolerance); 
    
            // Upstroke of cycle
            // -- Correction: u(h) <- u(h) + I(u(2h) - R.u(h))
            for (Level level=coarsestLevel+1; level<fineLevel; level++) {
                // Calculate I(u(2h) - R.u(h)), store in temporary 
                solution.coarsen(level, temp[level-1]);   // R.u(h))   
                solution[level-1] -= temp[level-1];       // u(2h) - R.u(h)
                solution.refine(level-1, temp[level]);    // I(u(2h) - R.u(h))    
                
                // Update u and do post-correction relaxation
                solution[level] += temp[level];
                relax(level, postRelax);
            }           
            
            // // Check whether residual has been reduced to truncation 
            // // error for the current cycle and break loop if so. 
            // evaluate_residual(fineLevel, temp[fineLevel]);             
            // if (temp[fineLevel].norm() < truncError) 
            //     break; 
        }                          
    }
}