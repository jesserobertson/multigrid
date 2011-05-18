/*
    fdvecarray.cpp (Multigrid)
    Jess Robertson, 2011-01-28
*/                            

#include "fdvecarray.hpp"

// Ctor   
mgrid::FDVecArray::FDVecArray(const double aspectRatio, const int nx, const int nz): 
    blitz::Array<double, 3>::Array(nx, nz, 2)
{
    calculate_geometry(aspectRatio, nx, nz);     
    (*this) = 0; // <-- Initialise data to 0
    
    // Generate component arrays             
    first.calculate_geometry(aspectRatio, nx, nz);
    second.calculate_geometry(aspectRatio, nx, nz);
    first.reference((*this)(blitz::Range::all(), blitz::Range::all(), 0)); 
    second.reference((*this)(blitz::Range::all(), blitz::Range::all(), 1));
}                                   

// Write method
void mgrid::FDVecArray::write(std::string fileString) {
    fileString.append(".nc");   // Add suffix to filename
    std::auto_ptr<NcFile> file(new NcFile(fileString.c_str(), NcFile::Replace));
    if (file->is_valid()) {
        // Define and add dimentions and variable
        NcDim* xDim = file->add_dim("x", nx);
        NcDim* zDim = file->add_dim("z", nz); 
        NcVar* xComp = file->add_var("u_x", ncDouble, xDim, zDim);
        NcVar* zComp = file->add_var("u_z", ncDouble, xDim, zDim);
        
        // Push data to file - data copied to ensure its contiguous and double
        // prescision, since NetCDF doesn't like anything bigger than a double
        blitz::Array<double, 2> xValues = (*this).first.copy();    
        blitz::Array<double, 2> zValues = (*this).second.copy();
        xComp->put(&xValues(0,0), nx, nz);        
        zComp->put(&zValues(0,0), nx, nz);
    }
}

// Magnitude, divergence etc...
void mgrid::FDVecArray::magnitude(ArrayType& result) {
    ARRAY_LOOP(result) 
        result(i, j) = sqrt(power<2>((*this)(i, j, 0)) + power<2>((*this)(i, j, 1)));                 
}  
void mgrid::FDVecArray::divergence(ArrayType& result) { 
    ARRAY_LOOP(result) 
        result(i, j) = first.dx(i, j) + second.dz(i, j);
}