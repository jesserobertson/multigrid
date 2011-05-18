/*
    fdvecarray.hpp (Multigrid)
    Jess Robertson, 2011-01-28
    
    Two-dimensional vector array class
*/

#ifndef FDVECARRAY_HPP_IF30D0T0
#define FDVECARRAY_HPP_IF30D0T0   
 
#include "fdbase.hpp"
#include "fdarray.hpp" 

namespace mgrid {

class FDVecArray: public FDBase, public blitz::Array<double, 3> {
public:
    FDVecArray(const double aspectRatio, const int nx, const int nz);  
    FDVecArray() {};
    virtual ~FDVecArray() {};        
    
    // Override base method for array resizing and boundary conditions
    inline virtual 
        void resize(const double aspectRatio, const int nx, const int nz); 
    inline virtual 
        void reference(FDVecArray array);
    inline virtual 
        void reference(blitz::Array<double, 3> array); 
    
    // Component arrays
    FDArray first, second;
    
    // Derivative methods
    void magnitude(ArrayType& result);
    void divergence(ArrayType& result);    
    
    // Writing method
    void write(std::string filestring);        
    
    // Explicitly inherit operators from blitz::Array, since 
    // these are not inherited by default.          
    using blitz::Array<double, 3>::operator=;
    using blitz::Array<double, 3>::operator+=;
    using blitz::Array<double, 3>::operator-=;
    using blitz::Array<double, 3>::operator*=;
    using blitz::Array<double, 3>::operator/=;      
    
    // Overload operators to use FDArray input. These employ the same templates
    // from blitz, but stop the compiler complaining about FDVecArray types.  
    inline FDVecArray& operator=(const FDVecArray& x);
    inline FDVecArray& operator+=(const FDVecArray& x);
    inline FDVecArray& operator-=(const FDVecArray& x);        
    inline FDVecArray& operator*=(const FDVecArray& x);
    inline FDVecArray& operator/=(const FDVecArray& x);
};

// Overloaded methods from blitz
inline void FDVecArray::resize(const double aspect, const int nx, const int nz) {
    calculate_geometry(aspect, nx, nz);
    blitz::Array<double, 3>::resize(nx, nz, 2);    
     
    // Generate component arrays                
    first.calculate_geometry(aspect, nx, nz);
    second.calculate_geometry(aspect, nx, nz);
    first.reference((*this)(blitz::Range::all(), blitz::Range::all(), 0)); 
    second.reference((*this)(blitz::Range::all(), blitz::Range::all(), 1));   
}            
inline void FDVecArray::reference(FDVecArray array) {
    calculate_geometry(array.aspectRatio, array.rows(), array.columns());
    blitz::Array<double, 3>::reference(array);
}
inline void FDVecArray::reference(blitz::Array<double, 3> array) {
    blitz::Array<double, 3>::reference(array);
}

// Overloaded operators for FDVecArray  
inline FDVecArray& FDVecArray::operator=(const FDVecArray& x) {
    using namespace blitz;
    (*this) = _bz_ArrayExpr<FastArrayIterator<double, 3> >(x.beginFast());
    return (*this);
}
inline FDVecArray& FDVecArray::operator+=(const FDVecArray& x) {
    using namespace blitz;
    (*this) += _bz_ArrayExpr<FastArrayIterator<double, 3> >(x.beginFast());
    return (*this);
}
inline FDVecArray& FDVecArray::operator-=(const FDVecArray& x) {
    using namespace blitz;
    (*this) -= _bz_ArrayExpr<FastArrayIterator<double, 3> >(x.beginFast());
    return (*this);
}        
inline FDVecArray& FDVecArray::operator*=(const FDVecArray& x) {
    using namespace blitz;
    (*this) *= _bz_ArrayExpr<FastArrayIterator<double, 3> >(x.beginFast());
    return (*this);
}
inline FDVecArray& FDVecArray::operator/=(const FDVecArray& x) {
    using namespace blitz;
    (*this) /= _bz_ArrayExpr<FastArrayIterator<double, 3> >(x.beginFast());
    return (*this);
}  

} // end namespace mgrid 

#endif /* end of include guard: FDVECARRAY_HPP_IF30D0T0 */
