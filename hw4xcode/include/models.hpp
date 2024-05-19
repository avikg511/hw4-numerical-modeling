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
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace ublas = boost::numeric::ublas;

class Models {
public:
    double dt;
    double dx;
    double duration;
    double rms;
    boost::numeric::ublas::matrix<double> data;
    Models(double timeRes, double dur, double spaceRes);
    void runForwardEuler();
    void runForwardEulerMatrix();
    void runImplicitEulerMatrix();
    void runMatsuno();
    void findRMSError();
    inline double exactSoln(int timeStep, int spaceStep);
    bool InvertMatrix(const ublas::matrix<double>& input, ublas::matrix<double>& inverse);
};

#endif /* models_hpp */
