/*
    fdarray.cpp (Multigrid)
    Jess Robertson, 2011-01-28     
    
    Two-dimensional scalar array class
*/             

#include "fdarray.hpp"
#include <netcdfcpp.h>

// Ctor
mgrid::FDArray::FDArray(const double aspectRatio, const int nx, const int nz):
    blitz::Array<double, 2>::Array(nx, nz) 
{
    calculate_geometry(aspectRatio, nx, nz); 
    *this = 0; // <-- Initialise data to 0
}                                       

// Boundary condition methods   
void mgrid::FDArray::update_boundaries() {
    // Loop over each boundary
    foreach(BoundaryFlag boundaryFlag, allBoundaryFlags) { 
        // Assign variable values depending on which boundary we are at 
        int dx, dz, sign; double spacing; blitz::Range i, j;  
        if (boundaryFlag == leftBoundary) {
            i       = blitz::Range(0); 
            j       = blitz::Range::all();
            sign    = -1; 
            dx      = 1; 
            dz      = 0; 
            spacing = hx;
        } else if (boundaryFlag == rightBoundary) {
            i       = blitz::Range(nx-1); 
            j       = blitz::Range::all();
            sign    = 1; 
            dx      = -1; 
            dz      = 0; 
            spacing = hx;
        } else if (boundaryFlag == topBoundary) {
            i       = blitz::Range::all(); 
            j       = blitz::Range(0);
            dx      = 0; 
            sign    = -1; 
            dz      = 1; 
            spacing = hz;
        } else if (boundaryFlag == bottomBoundary) {
            i       = blitz::Range::all(); 
            j       = blitz::Range(nz-1);
            dx      = 0; 
            sign    = 1; 
            dz      = -1; 
            spacing = hz;
        }
        
        // Actually perform update  
        foreach(BoundaryPoint pt, boundaryConditions.get(boundaryFlag))
            if (pt.conditionType == dirichlet) {
                (*this)(i, j) = pt.value;
            } else if (pt.conditionType == neumann) {
                (*this)(i, j) = (sign*12*(pt.value)*spacing 
                    + 48*(*this)(i+dx, j+dz) - 36*(*this)(i+2*dx, j+2*dz)
                    + 16*(*this)(i+3*dx, j+3*dz) - 3*(*this)(i+4*dx, j+4*dz)
                    )/25.0;
            } 
    }      
}

// Write method
void mgrid::FDArray::write(std::string fileString) {
    fileString.append(".nc");   // Add suffix to filename
    std::auto_ptr<NcFile> file(new NcFile(fileString.c_str(), NcFile::Replace));
    if (file->is_valid()) {
        // Define and add dimentions and variable
        NcDim* xDim = file->add_dim("x", nx);
        NcDim* zDim = file->add_dim("z", nz); 
        NcVar* variable = file->add_var("value", ncDouble, xDim, zDim);
        
        // Push data to file - data copied to ensure its contiguous and double
        // prescision, since NetCDF doesn't like anything bigger than a double
        blitz::Array<double, 2> values = (*this).copy();
        variable->put(&values(0,0), nx, nz);
    }
} 

// Gradient etc...
void mgrid::FDArray::gradient(VecArrayType& result) {
    ARRAY_LOOP((*this)) {
        result(i, j, 0) = (*this).dx(i, j);
        result(i, j, 1) = (*this).dz(i, j);
    }                  
}
void mgrid::FDArray::gradient_magnitude(ArrayType& result) {
    ARRAY_LOOP((*this)) 
        result(i, j) = sqrt(power<2>((*this).dx(i, j)) + power<2>((*this).dz(i, j)));
}

// Integrals
double mgrid::FDArray::calculate_flux() {  
    /*
        This uses an extended Simpsons rule, of the form:
        int(f(1..n)) = h*(17/48*f(0) + 59/48*f(1) + 43/48*f(2) + 49/48*f(3)
                + f(4) + f(5) + f(6) .... + f(n-5) + f(n-4)
                + 49/48*f(n-4) + 43/48*f(n-3) + 59/48*f(n-2) + 17/48*f(n-1))
    */
    // Interpolating polynomial coefficients (precalculated)
    static const double A = 17/48.;
    static const double B = 59/48.;
    static const double C = 43/48.;
    static const double D = 49/48.;  
      
    // Do integrals along array rows
    blitz::Array<double, 1> rowIntegrals(nx);
    for (int i=0; i<nx; i++) { 
        rowIntegrals(i) = A*(*this)(i, 0) + B*(*this)(i, 1) \
            + C*(*this)(i, 2) + D*(*this)(i, 3);              
        for (int j=4; j<nz-4; j++) rowIntegrals(i) += (*this)(i, j); 
        rowIntegrals(i) += A*(*this)(i, nz-1) + B*(*this)(i, nz-2) \
            + C*(*this)(i, nz-3) + D*(*this)(i, nz-4);   
        rowIntegrals(i) *= hz;
    }
    
    // Do integrals along array columns
    double result = A*rowIntegrals(0) + B*rowIntegrals(1) + C*rowIntegrals(2)
        + D*rowIntegrals(3);
    for (int i=4; i<nx-4; i++) result += rowIntegrals(i);
    result += A*rowIntegrals(nx-1) + B*rowIntegrals(nx-2) 
        + C*rowIntegrals(nx-3) + D*rowIntegrals(nx-4); 
    result *= hx;
    
    // Return result
    return result;
}
