//
//  models.hpp
//  hw4xcode
//
//  Created by Avik Ghosh on 5/11/24.
//

#ifndef models_hpp
#define models_hpp

#include <stdio.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/algorithm/cxx11/iota.hpp>
//#include <boost/config/warning_disable.hpp>

using namespace boost::numeric;

class Models {
public:
    double dt;
    double dx;
    double duration;
    double rms;
    boost::numeric::ublas::matrix<double> data;
    Models(double timeRes, double dur, double spaceRes);
    void runForwardEuler();
    void runMatsuno();
    void findRMSError();
    inline double exactSoln(int timeStep, int spaceStep);
};

#endif /* models_hpp */
