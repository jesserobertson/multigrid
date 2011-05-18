/*
    multigrid_expections (Multigrid)
    Jess Robertson, 2011-01-28
    
    Basic exception classes for multigrid
*/                                       

#ifndef MULTIGRID_EXCEPTIONS_HPP_6MQVVDHT
#define MULTIGRID_EXCEPTIONS_HPP_6MQVVDHT

#include <iostream>

namespace mgrid { 
    
enum MessageType {ErrorMessage, WarningMessage, StatusMessage};

class Message: public std::ostringstream {
public:
    Message(MessageType messageType) {
        switch (messageType) {
            case ErrorMessage:  
                (*this) << " *ERROR* ";
                break;
            case WarningMessage:
                (*this) << " Warning ";  
                break;
            case StatusMessage:   
                (*this) << "  ---->  ";  
                break;
        }        
    }
};  

class MultigridException: public std::exception {};

class BoundaryLengthException: public MultigridException { 
public:
    virtual const char* what() const throw() {
        Message msg(ErrorMessage); 
        msg << "Lengths of arrays supplied to boundary conditions "
            << "are different.";
        return msg.str().c_str();
    }
};  

class InvalidBoundaryCondition: public MultigridException {
public:
    virtual const char* what() const throw() {
        Message msg(ErrorMessage); 
        msg << "A boundary condition of the wrong length was supplied.";
        return msg.str().c_str();
    }
};  

} // end namespace mgrid



#endif /* end of include guard: MULTIGRID_EXCEPTIONS_HPP_6MQVVDHT */
