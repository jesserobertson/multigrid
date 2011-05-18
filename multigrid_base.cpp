/*
    multigrid_base.cpp (Multigrid)
    Jess Robertson, 2011-01-28 
    
    Base class for multigrid solvers
*/                            

#include "multigrid_base.hpp"

// Ctor
mgrid::MultigridBase::MultigridBase(const Settings& settings):  
    solution(settings),
    source(settings),
    temp(settings),       
    cycleType(settings.mgCycleType),
    preRelax(settings.preMGRelaxIter), 
    postRelax(settings.postMGRelaxIter), 
    residualTolerance(settings.residualTolerance), 
    maxIterations(settings.maximumIterations),  
    aspect(settings.aspectRatio),
    sourceIsSet(false),
    initialIsSet(false)
{
    // Initialise some other variables
    finestLevel = solution.finestLevel;
    coarsestLevel = solution.coarsestLevel;  
    nxfine = solution[finestLevel].rows(); 
    nzfine = solution[finestLevel].columns(); 
}

// Relaxation methods
void mgrid::MultigridBase::relax(const Level level, const unsigned long N) {
    // Relax for N iterations
    for (unsigned long iter=0; iter<N; iter++) { 
        RED_BLACK_LOOP(solution[level])
            relaxation_updater(level, i, j); 
        
        // Update boundaries
        solution[level].update_boundaries();
    }      
}                                       
void mgrid::MultigridBase::relax(const Level level, const double tolerance) {
    // Relax until specified tolerance
    for (unsigned long iter=0; iter<maxIterations; iter++) {  
        double residualSum = 0, normSum = 0, tmp;
        RED_BLACK_LOOP(solution[level]) {
     	    // Store current value, calculate update  
            tmp = solution[level](i, j);    
            relaxation_updater(level, i, j);    
                                                       
            // Calculate change and add to sum
            tmp = (solution[level](i, j) - tmp);          
            residualSum += power<2>(tmp);
            normSum += power<2>(solution[level](i, j));  
        }            
     	
     	// Update boundaries                             
        solution[level].update_boundaries();                                                           
     	
     	// Check for convergence
     	if ((sqrt(residualSum)/sqrt(normSum)) < tolerance) return; 
    } 
}                      

// Write method
void mgrid::MultigridBase::write(int numOfVariables, std::string fileRoot) { 
    // Get generated file name from settings instance
    fileRoot.append(filename()).append(".nc");      
    
    // Generate new netCDF file using given fileString - replace the file if
    // it already exists
    std::auto_ptr<NcFile> file(new NcFile(fileRoot.c_str(), NcFile::Replace));
    
    // Check to see that file has been generated  
    if (file->is_valid()) {
        // Define and add dimentions and variables
        NcDim* xDim = file->add_dim("x", nxfine);
        NcDim* zDim = file->add_dim("z", nzfine);              
        NcVar* XxVals = file->add_var("x", ncDouble, xDim);   
        NcVar* ZzVals = file->add_var("z", ncDouble, zDim);
        const double hx = solution[finestLevel].spacing(0);
        const double hz = solution[finestLevel].spacing(1);
        double XxArray[nxfine];
        double ZzArray[nzfine];
        for (int i=0; i<nxfine; i++) XxArray[i] = i*hx;
        for (int j=0; j<nzfine; j++) ZzArray[j] = j*hz; 
        XxVals->put(XxArray, nxfine);
        ZzVals->put(ZzArray, nzfine);
        
        // Add solution
        NcVar* variable = file->add_var("solution", ncDouble, xDim, zDim);
        variable->put(&solution[finestLevel](0,0), nxfine, nzfine); 

        // Calculate gradient components
        if (numOfVariables > 1) {
            NcVar* modGrad = file->add_var("gradient", ncDouble, xDim, zDim);  
            FDArray gradientValues(aspect, nxfine, nzfine); 
            solution[finestLevel].gradient_magnitude(gradientValues);     
            modGrad->put(&gradientValues(0,0), nxfine, nzfine);
        } 
        
        // Calculate residual values
        if (numOfVariables > 2) {
            NcVar* resid = file->add_var("log_residual", ncDouble, xDim, zDim);
            FDArray residualValues(aspect, nxfine, nzfine);
            evaluate_residual(finestLevel, residualValues); 
            residualValues = log10(residualValues);       
            resid->put(&residualValues(0,0), nxfine, nzfine); 
        } 
        
        // Add some attributes describing settings 
        file->add_att("aspect_ratio", aspect);
    }
}