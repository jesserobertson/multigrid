/*
    bingham.cpp (Multigrid)
    Jess Robertson, 2010-08-09
*/

#include "mosolov.hpp"

using namespace mgrid;

Mosolov::Mosolov(const MosolovSettings& settings):
    LinearMultigrid::LinearMultigrid(settings.multigridSettings),
    temp(settings.multigridSettings.aspectRatio, nxfine, nzfine), 
    multiplier(settings.multigridSettings.aspectRatio, nxfine, nzfine),
    strainRate(settings.multigridSettings.aspectRatio, nxfine, nzfine),
    determinant(settings.multigridSettings.aspectRatio, nxfine, nzfine),
    tempVec(settings.multigridSettings.aspectRatio, nxfine, nzfine),
    solutionGradient(settings.multigridSettings.aspectRatio, nxfine, nzfine),  
    aspectRatio(settings.multigridSettings.aspectRatio),
    binghamNumber(settings.binghamNumber),
    maxLagrangeIteration(settings.maxLagrangeIteration),
    lagrangeTolerance(settings.lagrangeTolerance),
    alpha(settings.augmentingParameter)
{
    // Set boundary conditions for velocity array 
    solution.boundaryConditions.set(leftBoundary,   zeroNeumannCondition);
    solution.boundaryConditions.set(rightBoundary,  zeroDirichletCondition);
    solution.boundaryConditions.set(topBoundary,    zeroNeumannCondition);    
    solution.boundaryConditions.set(bottomBoundary, zeroDirichletCondition); 
    
    // Set source and multiplier values   
    source[finestLevel] = -1.0; 
    sourceIsSet = true;     
    multiplier = 0.0;  
    
    // Specify that problem has been set up on given thread
    std::ostringstream msg;
    msg << " -- Problem (" << aspectRatio << ", " << binghamNumber 
        << ") initialised." << std::endl;
    std::cout << msg.str(); std::cout.flush();
}    
Mosolov::~Mosolov() { /* pass */ }

void Mosolov::write(int numOfVariables, std::string fileRoot) {
    // Get generated file name from settings instance
    fileRoot.append(filename()).append(".nc");      
    
    // Generate new netCDF file using given fileString - replace the file if
    // it already exists. File will sync and close automatically when it 
    // goes out of scope.
    std::auto_ptr<NcFile> file(new NcFile(fileRoot.c_str(), NcFile::Replace));
    
    // Check to see that file has been generated
    if (file->is_valid()) {           
        // Define and add attributes which specify aspect ratio, Bingham 
        // number and total flux 
        file->add_att("bingham_number", binghamNumber);
        file->add_att("aspect_ratio", aspectRatio);
        file->add_att("total_flux", solution[solution.finestLevel].calculate_flux());
        
        // Specify axis values   
        const double hx = solution[solution.finestLevel].spacing(0);
        const double hz = solution[solution.finestLevel].spacing(1);  
        double XxArray[nxfine];
        double ZzArray[nzfine];
        for (int i=0; i<nxfine; i++) XxArray[i] = i*hx;
        for (int j=0; j<nzfine; j++) ZzArray[j] = j*hz;
        
        // Define and add dimentions and variables
        NcDim* Xx = file->add_dim("x", nxfine);
        NcDim* Zz = file->add_dim("z", nzfine);  
        NcVar* XxVals = file->add_var("x", ncDouble, Xx);   
        NcVar* ZzVals = file->add_var("z", ncDouble, Zz);
        NcVar* variable = file->add_var("velocity", ncDouble, Xx, Zz);
        
        // Add axis values and velocity components 
        XxVals->put(XxArray, nxfine);
        ZzVals->put(ZzArray, nzfine);
        variable->put(&solution[solution.finestLevel](0,0), nxfine, nzfine);   
        
        // Calculate and add gradient components
        if (numOfVariables > 1) {    
            NcVar* modGrad = file->add_var("strain_rate", ncDouble, Xx, Zz); 
            solution[solution.finestLevel].gradient_magnitude(temp);     
            modGrad->put(&temp(0,0), nxfine, nzfine); 
        }
        
        // Calculate and add residual
        if (numOfVariables > 2) {    
            NcVar* resid = file->add_var("log_residual", ncDouble, Xx, Zz);
            evaluate_residual(finestLevel, temp); 
            temp = log10(temp);       
            resid->put(&temp(0,0), nxfine, nzfine);   
        }   
        
        // Let std::cout know that the file has been written
        std::ostringstream msg;
        msg << " -- Problem (" << aspectRatio << ", " << binghamNumber 
            << ") solution written to " << fileRoot << std::endl;
        std::cout << msg.str(); std::cout.flush(); 
    }
}

void Mosolov::solve() {
    // Solve initial problem, then loop through augmented Lagrangian iteration
    multigrid();    
    for(unsigned int iter = 0; iter < maxLagrangeIteration; ++iter) {
        // Calculate new strain rate  
        solution[finestLevel].gradient(solutionGradient);
        determinant = alpha*solutionGradient + multiplier;
        ARRAY_LOOP(strainRate.first) {           
            const double d1sq = determinant.first(i, j)*determinant.first(i, j);
            const double d2sq = determinant.second(i, j)*determinant.second(i, j);
            const double detMagnitude = sqrt(d1sq + d2sq);
            if (detMagnitude*detMagnitude <= binghamNumber*binghamNumber) {
                strainRate.first(i, j) = 0; 
                strainRate.second(i, j) = 0;
            } else {                       
                strainRate.first(i, j) = (1-binghamNumber/detMagnitude)
                    *determinant.first(i, j)/alpha;
                strainRate.second(i, j) = (1-binghamNumber/detMagnitude)
                    *determinant.second(i, j)/alpha;
            }   
        } 
        
        // Construct new right hand side
        tempVec.first = alpha*strainRate.first - multiplier.first;
        tempVec.second = alpha*strainRate.second - multiplier.second;
        tempVec.divergence(temp);
        source[finestLevel] = (temp - 1.0)/(1+alpha);
        
        // Calculate new velocity
        multigrid();        
        
        // Check for convergence (i.e. when $\dot\gamma = \nabla u$)
        const double resid = _normed_residual();   
        if (resid < lagrangeTolerance) {   
            std::ostringstream msg;
            msg << " -- Problem (" << aspectRatio << ", " << binghamNumber 
                << ") converged after " << iter << " iterations. " << std::endl 
                <<"    Residual = " << resid << std::endl; 
            std::cout << msg.str(); std::cout.flush();
            return;
        }
        
        // Calculate new multiplier                  
        solution[solution.finestLevel].gradient(solutionGradient);
        multiplier.first += 
            alpha*(solutionGradient.first - strainRate.first);
        multiplier.second +=
            alpha*(solutionGradient.second - strainRate.second);
    } 
    
    // If we're here, then the convergence has failed  
    std::ostringstream msg;
    msg << " -- Problem (" << aspectRatio << ", " << binghamNumber 
        << ") failed to converge after " << maxLagrangeIteration
        << " iterations. Residual = " << _normed_residual() << std::endl;
    std::cout << msg.str(); std::cout.flush();
}  