/*  
    utilities.hpp (Multigrid)
    Jess Robertson, 2010-07-25   
    
    Some useful templates for multigrid solver
*/

#ifndef TEMPLATES_HPP_PM2L3HHD
#define TEMPLATES_HPP_PM2L3HHD

#include <string>   
#include <fstream>
#include <iostream>
#include <math.h>
#include <blitz/array.h>
#include <boost/tuple/tuple.hpp>  
#include <boost/foreach.hpp>     

#include "types.hpp"      

namespace mgrid {
    
// Ensure that blitz++ uses threadsafe array access
#define BZ_THREADSAFE   
    
// Redefine boost_foreach to look a bit nicer
#ifndef foreach
#define foreach BOOST_FOREACH
#endif

// Defines a loop to loop over the given array
#ifndef ARRAY_LOOP
#define ARRAY_LOOP(array) \
    for (int i=0; i < array.rows(); i++) \
        for (int j=0; j < array.columns(); j++)
#endif    

// Defines a loop which loops over red points and then black points   
#ifndef RED_BLACK_LOOP
#define RED_BLACK_LOOP(array) \
    int jof = 1, iof = 1; \
    for (int pColor=1; pColor<3; pColor++, iof=3-jof, jof=3-jof) \
        for (int j=1; j < array.columns()-1; j++, iof=3-iof) \
         	for (int i=iof; i < array.rows()-1; i+=2)
#endif                    

// Defines a dot product for two vector arrays
#ifndef dot_product
#define dot_product(A, B) (A.first*B.first + A.second*B.second)
#endif    

// Template for to expand integer power expressions
template <int n> inline double power(const double& m) { 
    return power<n-1>(m)*m; 
}         
template <> inline double power<1>(const double& m) { return m; }  // x^1 = x
template <> inline double power<0>(const double& m) { return 1; }  // x^0 = 1

// Template to turn any class with a << operator into a string
template <class T> std::string str(T i) {
    std::ostringstream buffer; 
	buffer << i;
    return buffer.str();
}        

// Bool printer for logging messages
inline std::string bool_string(const bool& arg) {
	if (arg) return "true";
	else return "false";
}             

// Template for transferring two-member tiny vector to a tuple and once for 
// making a tinyvector
template <typename T>
inline boost::tuple<T, T> to_tuple(blitz::TinyVector<T, 2> input) {
    return boost::make_tuple(input(0), input(1));
}                                 
template <typename T>
inline blitz::TinyVector<T, 2> to_tinyvec(T first, T second) {
    blitz::TinyVector<T, 2> returnVec;
    returnVec(0) = first;
    returnVec(1) = second;
    return returnVec;
}     

} // end namespace mgrid  

#endif /* end of include guard: TEMPLATES_HPP_PM2L3HHD */