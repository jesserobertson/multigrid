/*
    FDArray.hpp (Multigrid)
    Jess Robertson, 2011-01-28     
    
    Two-dimensional scalar array class
*/                            

#ifndef FDARRAY_HPP_2W1G76VE
#define FDARRAY_HPP_2W1G76VE

#include <boost/foreach.hpp>     
#include <boost/function.hpp>      
#include <boost/bind.hpp>  

#include "fdbase.hpp"  
#include "boundary_conditions.hpp"

namespace mgrid {

class FDArray: public FDBase, public blitz::Array<double, 2> {
public:
    FDArray(const double aspectRatio, const int nx, const int nz);  
    FDArray() {}; 
    virtual ~FDArray() {};             

    // Boundary condition methods
    BoundaryConditions boundaryConditions;
    void update_boundaries();

    // Writing method
    virtual void write(std::string filestring);            

    // Override base methods for array resizing and referencing
    inline virtual 
        void resize(const double aspectRatio, const int nx, const int nz);    
    inline virtual 
        void reference(blitz::Array<double, 2> array);     
    inline virtual 
        void reference(FDArray array);       

    // Gradient etc...
    void gradient(VecArrayType& result);       
    void gradient_magnitude(ArrayType& result);   

    // Integration
    double calculate_flux();  

    // Add a simple function to calculate the Frobenius norm of an array
    inline double norm();

    // Explicitly inherit operators from blitz::Array, since 
    // these are not inherited by default.          
    using blitz::Array<double, 2>::operator=;
    using blitz::Array<double, 2>::operator+=; 
    using blitz::Array<double, 2>::operator-=;
    using blitz::Array<double, 2>::operator*=;
    using blitz::Array<double, 2>::operator/=;      

    // Overload operators to use FDArray input. These employ the same templates
    // from blitz, but stop the compiler complaining about FDArray types.  
    inline FDArray& operator=(const FDArray& x); 
    inline FDArray& operator=(const blitz::Array<double, 2>& x);
    inline FDArray& operator+=(const FDArray& x);
    inline FDArray& operator-=(const FDArray& x);        
    inline FDArray& operator*=(const FDArray& x);
    inline FDArray& operator/=(const FDArray& x); 

    // Settor functions which take a function as input
    inline FDArray& operator=(const boost::function<double (double, double)> f);

    // Derivatives methods
    inline void derivatives(Deriv& du, const int i, const int j); 
    inline double dx(const int i, const int j);
    inline double dxu(const int i, const int j);
    inline double dz(const int i, const int j);      
    inline double dzu(const int i, const int j);    
    inline double dxx(const int i, const int j);
    inline double dxxu(const int i, const int j);
    inline double dzz(const int i, const int j);
    inline double dzzu(const int i, const int j);
    inline double dxz(const int i, const int j);
    inline double dxzu(const int i, const int j);  

    // Array-wide derivative methods
    inline void dx(ArrayType& result);
    inline void dz(ArrayType& result);
    inline void dxx(ArrayType& result);
    inline void dzz(ArrayType& result);    
    inline void dxz(ArrayType& result);             
};                

// Overloaded methods from Array
inline void FDArray::resize(const double aspectRatio, const int nx, const int nz) {
    calculate_geometry(aspectRatio, nx, nz);
    blitz::Array<double, 2>::resize(nx, nz);
    boundaryConditions.resize(nx, nz);        
}      
inline void FDArray::reference(FDArray referent) {
    calculate_geometry(referent.aspectRatio, referent.rows(), referent.columns());
    boundaryConditions.reference(referent.boundaryConditions);
    blitz::Array<double, 2>::reference(referent);
}
inline void FDArray::reference(blitz::Array<double, 2> referent) {
    blitz::Array<double, 2>::reference(referent);
}      

// Overloaded array operators for FDArray
inline FDArray& FDArray::operator=(const FDArray& x) {
    using namespace blitz;
    (*this) = _bz_ArrayExpr<FastArrayIterator<double, 2> >(x.beginFast());
    return (*this);
}  
inline FDArray& FDArray::operator=(const blitz::Array<double, 2>& x) {
    using namespace blitz;
    (*this) = _bz_ArrayExpr<FastArrayIterator<double, 2> >(x.beginFast());
    return (*this);
}
inline FDArray& FDArray::operator=(boost::function<double (double, double)> f) {
    ARRAY_LOOP((*this))
        (*this)(i, j) = f(i*hx, j*hz);   
    return (*this);  
}
inline FDArray& FDArray::operator+=(const FDArray& x) {
    using namespace blitz;
    (*this) += _bz_ArrayExpr<FastArrayIterator<double, 2> >(x.beginFast());
    return (*this);
}
inline FDArray& FDArray::operator-=(const FDArray& x) {
    using namespace blitz;
    (*this) -= _bz_ArrayExpr<FastArrayIterator<double, 2> >(x.beginFast());
    return (*this);
}        
inline FDArray& FDArray::operator*=(const FDArray& x) {
    using namespace blitz;
    (*this) *= _bz_ArrayExpr<FastArrayIterator<double, 2> >(x.beginFast());
    return (*this);
}
inline FDArray& FDArray::operator/=(const FDArray& x) {
    using namespace blitz;
    (*this) /= _bz_ArrayExpr<FastArrayIterator<double, 2> >(x.beginFast());
    return (*this);
}     

// FDArray derivatives calculation methods 
inline void FDArray::derivatives(Deriv& du, const int i, const int j) {
    du.dx = (*this).dx(i, j);
    du.dxu = (*this).dxu(i, j);
    du.dz = (*this).dz(i, j);
    du.dzu = (*this).dzu(i, j);
    du.dxx = (*this).dxx(i, j);
    du.dxxu = (*this).dxxu(i, j);
    du.dzz = (*this).dzz(i, j);
    du.dzzu = (*this).dzzu(i, j);
    du.dxz = (*this).dxz(i, j);
    du.dxzu = (*this).dxzu(i, j);
}
inline double FDArray::dx(const int i, const int j) {
    if (i == 0) {
		// forward difference in i
        return (3*(*this)(i,j) - 4*(*this)(i+1,j) + (*this)(i+2,j))*xfactor;
	} else if (i == nx-1) {
	    // backward difference in i     
		return (-(*this)(i-2,j) + 4*(*this)(i-1,j) - 3*(*this)(i,j))*xfactor;
	} else {
	    // centered difference in i
		return ((*this)(i-1,j) - (*this)(i+1,j))*xfactor;
	}
} 
inline double FDArray::dxu(const int i, const int j) {
    if (i == 0) { 
		// forward difference in i
	    return 3*xfactor;
	} else if (i == nx-1) {
	    // backward difference in i     
		return -3*xfactor;
	} else {
	    // centered difference in i
		return 0;
	}
}
inline double FDArray::dz(const int i, const int j) {
    if (j == 0) { 
	    // forward difference in j 
		return (-4*(*this)(i,j+1) + 3*(*this)(i,j) + (*this)(i,j+2))*zfactor;
	} else if (j == nz-1) {
	    // backward difference in j 
		return (-(*this)(i,j-2) + 4*(*this)(i,j-1) - 3*(*this)(i,j))*zfactor;
	} else {
	    // centered difference in j
		return -((*this)(i,j+1) - (*this)(i,j-1))*zfactor;
	}
}
inline double FDArray::dzu(const int i, const int j) {
    if (j == 0) { 
	    // forward difference in j 
		return 3*zfactor;
	} else if (j == nz-1) {
	    // backward difference in j 
        return -3*zfactor;	    
	} else {
	    // centered difference in j
		return 0;
	}

}
inline double FDArray::dxx(const int i, const int j) {
    if (i == 0) { 
		// forward difference in i   
		return (-(*this)(i+3,j) + 4*(*this)(i+2,j) - 5*(*this)(i+1,j) 
		    + 2*(*this)(i,j))*xxfactor;
	} else if (i == nx-1) {         
		// backward difference in i     
		return (-(*this)(i-3, j) + 4*(*this)(i-2,j) - 5*(*this)(i-1,j) 
		    + 2*(*this)(i,j))*xxfactor;
 	} else {
		// centered difference in i 
		return ((*this)(i-1,j) - 2*(*this)(i,j) + (*this)(i+1,j))*xxfactor;
	}	
}
inline double FDArray::dxxu(const int i, const int j) {
    if (i == 0) { // Calculate derivatives in the x direction 
		// forward difference in i   
		return 2*xxfactor;
	} else if (i == nx-1) {         
		// backward difference in i     
		return 2*xxfactor;
 	} else {
		// centered difference in i 
		return -2*xxfactor;
	}
}
inline double FDArray::dzz(const int i, const int j) {
    if (j == 0) { 
	    // forward difference in j     
		return (-(*this)(i,j+3)+ 4*(*this)(i,j+2) - 5*(*this)(i,j+1) 
		    + 2*(*this)(i,j))*zzfactor;
	} else if (j == nz-1) {
	    // backward difference in j  
		return (-(*this)(i, j-3) + 4*(*this)(i,j-2) - 5*(*this)(i,j-1) 
		    + 2*(*this)(i,j))*zzfactor;
	} else {    
	    // centered difference in j    
		return ((*this)(i,j-1) - 2*(*this)(i,j) + (*this)(i,j+1))*zzfactor;
	}       
}
inline double FDArray::dzzu(const int i, const int j) {
    if (j == 0) { 
	    // forward difference in j     
		return 2*zzfactor;
	} else if (j == nz-1) {
	    // backward difference in j  
		return 2*zzfactor;
	} else {    
	    // centered difference in j    
		return -2*zzfactor;
	}
}
inline double FDArray::dxz(const int i, const int j) {
    if (j == 0) {
		if (i == 0) {
			// forward difference in i, forward difference in j
			return (16*(*this)(i+1,j+1) - 12*((*this)(i+1,j) + (*this)(i,j+1)) 
				+ 9*(*this)(i,j) - 4*((*this)(i+1,j+2) + (*this)(i+2,j+1))
				+ 3*((*this)(i+2,j) + (*this)(i,j+2)) + (*this)(i+2,j+2))*xzfactor;    
		} else if (i == nx-1) {
			// backward difference in i, forward difference in j
			return (-16*(*this)(i-1,j+1) + 12*((*this)(i-1,j) + (*this)(i,j+1)) 
				- 9*(*this)(i,j) + 4*((*this)(i-1,j+2) + (*this)(i-2,j+1)) 
				- 3*((*this)(i-2,j) - (*this)(i,j+2)) - (*this)(i-2,j+2))*xzfactor;
		} else {
			// centered difference in i, forward difference in j
			return (4*((*this)(i+1,j+1)-(*this)(i-1,j+1)) 
				+ 3*((*this)(i-1,j)-(*this)(i+1,j)) + (*this)(i-1,j+2) 
				- (*this)(i+1,j+2))*xzfactor;      
		}
	} else if (j == nz-1) {
		if (i == 0) {
			// forward difference in i, backward difference in j
            return (-16*(*this)(i+1,j-1) + 12*((*this)(i+1,j) + (*this)(i,j-1)) 
				- 9*(*this)(i,j) + 4*((*this)(i+2,j-1) + (*this)(i+1,j-2)) 
				- 3*((*this)(i+2,j) + (*this)(i,j-2)) - (*this)(i+2,j-2))*xzfactor;  
		} else if (i == nx-1) {
			// backward difference in i, backward difference in j
			return (16*(*this)(i-1,j-1) - 12*((*this)(i-1,j) + (*this)(i,j-1)) 
				+ 9*(*this)(i,j) - 4*((*this)(i-1,j-2) + (*this)(i-2,j-1)) 
				+ 3*((*this)(i-2,j) + (*this)(i,j-2)) + (*this)(i-2,j-2))*xzfactor;      
		} else {
			// centered difference in i, backward difference in j
			return (4*((*this)(i-1,j-1) - (*this)(i+1,j-1)) 
			    + 3*((*this)(i+1,j) - (*this)(i-1,j)) 
			    - (*this)(i-1,j-2) + (*this)(i+1,j-2))*xzfactor;    
		}
	} else {
		if (i == 0) {
			// forward difference in i, centered difference in j
            return (3*((*this)(i,j-1) - (*this)(i,j+1))
				- 4*((*this)(i+1,j-1) - (*this)(i+1,j+1)) 
				+ (*this)(i+2,j-1) - (*this)(i+2,j+1))*xzfactor;   
		} else if (i == nx-1) {
			// backward difference in i, centered difference in j
			return (3*((*this)(i,j+1) - (*this)(i,j-1)) 
    			- 4*((*this)(i-1,j+1) - (*this)(i-1,j-1))
				+ (*this)(i-2,j+1) - (*this)(i-2,j-1))*xzfactor;         
		} else {
			// centered difference in i, centered difference in j
			return ((*this)(i-1,j-1) - (*this)(i-1,j+1) - (*this)(i+1,j-1) 
				+ (*this)(i+1,j+1))*xzfactor;      
		}
	}
}
inline double FDArray::dxzu(const int i, const int j) {
    if (j == 0) {
		if (i == 0) {
			// forward difference in i, forward difference in j   
			return 9*xzfactor;
		} else if (i == nx-1) {
			// backward difference in i, forward difference in j
			return -9*xzfactor;
		} else {
			// centered difference in i, forward difference in j      
			return 0;
		}
	} else if (j == nz-1) {
		if (i == 0) {
			// forward difference in i, backward difference in j 
			return -9*xzfactor;
		} else if (i == nx-1) {
			// backward difference in i, backward difference in j      
			return 9*xzfactor;
		} else {
			// centered difference in i, backward difference in j 
			return 0;
		}
	} else {
		if (i == 0) {
			// forward difference in i, centered difference in j  
			return 0;
		} else if (i == nx-1) {
			// backward difference in i, centered difference in j    
			return 0;
		} else {
			// centered difference in i, centered difference in j     
			return 0;
		}
	}
}

// Array-wide derivatives methods
inline void FDArray::dx(ArrayType& result) {
    ARRAY_LOOP((*this)) 
        result(i, j) = (*this).dx(i, j);
}
inline void FDArray::dz(ArrayType& result) {
    ARRAY_LOOP((*this)) 
        result(i, j) = (*this).dz(i, j);
}
inline void FDArray::dxx(ArrayType& result) {
    ARRAY_LOOP((*this)) 
        result(i, j) = (*this).dxx(i, j);
}
inline void FDArray::dzz(ArrayType& result) {
    ARRAY_LOOP((*this)) 
        result(i, j) = (*this).dzz(i, j);
}
inline void FDArray::dxz(ArrayType& result) {
    ARRAY_LOOP((*this)) 
        result(i, j) = (*this).dxz(i, j);
}

// Other inline functions
inline double FDArray::norm() {
    return sqrt(blitz::sum((*this)*(*this))/blitz::dot(shape(), shape()));
}       

} // end namespace mgrid

#endif /* end of include guard: FDARRAY_HPP_2W1G76VE */
