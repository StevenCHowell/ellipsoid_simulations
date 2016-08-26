/*
 * Modified and adpated into SASSIE by Hailiang Zhang in April 2014
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
 * Please note, that Artistic License is slightly more restrictive
 * then GPL license in terms of distributing the modified versions
 * of this software (they should be approved first).
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 * for details. Such license fits scientific software better then
 * GPL because it prevents the distribution of bugged derivatives.
 *
*/

#ifndef SASPTEROS_ERROR_H
#define SASPTEROS_ERROR_H

#include <sstream>
#include <iostream>
#include <exception>

using namespace std;

class SasPteros_error {
public:

    SasPteros_error(const SasPteros_error& p){
        text.str(p.text.str());
    }

    SasPteros_error(){
        text.str("");
    }

    /// Constructs an exception object with text message
    SasPteros_error(string s){
        text.str(s);
    }

    /** Operator << allows constructing error strings on the fly
         like this:
            \code throw SasPteros_error("Wrong number ") << 5
                    << " should be between "
                    << 7 << " and " << 10;
            \endcode
        */
    template<class T>
    SasPteros_error& operator<<(T data){
        text << data; //Collect data
        return *this;
    };

    /// Print error message
    void print(){
        cout << endl << "PTEROS terminated due to the following error:"
             << endl << text.str() << endl;
    }

    /// Return error message as string
    std::string what(){
        return text.str();
    }

private:
    std::stringstream text;
};

#endif
