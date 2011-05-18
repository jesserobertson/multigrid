/*
    multigrid_linear.cpp (Multigrid)
    Jess Robertson, 2011-01-28
*/                        

#include "multigrid_linear.hpp"

void mgrid::LinearMultigrid::multigrid() {
    // Check that source array has been set  
    if (not(sourceIsSet)) return;
    
    // Constants
    const Level finestLevel = solution.finestLevel;
    const Level coarsestLevel = solution.coarsestLevel;
    
    // Initialise right-hand-side
    for (Level level=finestLevel; level>0; level--) {
        source.coarsen(level);     
        solution.coarsen(level);
    }

    // Solve on coarsest level
    relax(coarsestLevel, residualTolerance);

    // Full Multigrid loop
    for (Level fineLevel=1; fineLevel<=finestLevel; fineLevel++) {
        // V-cycle loop at each (successively finer) level
        solution.refine(fineLevel-1); // interpolate to next level    
        for (int cycle=0; cycle < cycleType; cycle++) {
            // Downstroke of cycle:    
            //  -- New residual: r(2h) = 0 (see note below)
            //  -- New rhs: f(2h) = R.f(h) - L.u(2h)    
            // Note that each sucessive level in the downstroke is calculating
            // a _residual_, not a coarser version of the solution. Each level
            // therefore needs to be set to zero on the way down.
            for (Level level=fineLevel; level>0; level--) {
                relax(level, preRelax);
                evaluate_residual(level, temp[level]); 
                temp.coarsen(level, source[level-1]);     
                solution[level-1] = 0; // initialise next level's residual
            }

            // Solve problem on coarsest level    
            relax(coarsestLevel, residualTolerance);

            // Upstroke of cycle:
            // -- Correction: u(h) <- u(h) + I.u(2h)
            for (Level level=1; level<=fineLevel; level++) {
                solution.refine(level-1, temp[level]);     
                solution[level] += temp[level];
                relax(level, postRelax); 
            } 
        }
    } 
    
    // Do final update 
    relax(finestLevel, postRelax);
    solution[finestLevel].update_boundaries();
}
