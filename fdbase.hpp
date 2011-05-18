/*
    fdbase.hpp (Multigrid)
    Jess Robertson, 2011-01-28
*/    

#ifndef FDBASE_HPP_EK1M7AR9
#define FDBASE_HPP_EK1M7AR9

#include <netcdfcpp.h>

#include "types.hpp" 
#include "multigrid_exceptions.hpp"
#include "utilities.hpp"

namespace mgrid { 
                                                                                    
// = FDBase - abstract array class =     
class FDBase {
public:
    FDBase();
    virtual ~FDBase() {};    
    
    // Function to calculate spacing and resolution
    void calculate_geometry(const double aspectRatio, const int nx, const int nz);
    
    // Override base methods for array resizing
    virtual void resize(double aspectRatio, const int nx, const int nz)=0; 
        
    // Accessors
    inline const double spacing(const int n) {
        if (n == 0) return hx;
        else return hz; 
    }   
                    
protected:
    // Geometry attributes
    int nx, nz;      
    double aspectRatio, hx, hz;
    Deriv du;         
    double xfactor, zfactor;              // Denominators (like 1/hx) used  
	double xxfactor, zzfactor, xzfactor;  // for calculating derivatives  
};       

} // end namespace mgrid

#endif /* end of include guard: FDBASE_HPP_EK1M7AR9 */
