//
//  models.cpp
//  hw4xcode
//
//  Created by Avik Ghosh on 5/11/24.
//

#include "models.hpp"
Models::Models(double timeRes, double dur, double spaceRes) {
    
    this->dt = timeRes;
    this->dx = spaceRes;
    this->duration = dur;
    
    this->data = ublas::matrix<double>(int(ceil(this->duration / timeRes)), int(ceil(1 / spaceRes)));
    
//    Initial Condition: C(t=0, x) = 1 + cos (2 \pi * (dx * col)), where dx * col = x.
    for (int col = 0; col < data.size2(); col++) {
        data(0, col) = 1 + cos(2 * M_PI * (this->dx * col));
    }
}


void Models::runForwardEuler() {
    /*
     Forward Euler Code:
     \frac{\partial C}{ \partial t} = \frac{ \partial^2 C }{\partial x^2} \\
                        = \frac{1}{\Delta x^2} (C^n_{k+1} -2C^n_{k} + C^n_{k-1}) \\
     so C^{n+1}_k = C^n_{k} + \frac{ \Delta t }{\Delta x^2} (C^n_{k+1} -2C^n_{k} + C^n_{k-1})
     
     Initial Conditions: C(t = 0, x) = 1 + cos(2 \pi x)
     
     Exact Solution: C(t,x) = 1 + exp (-4 \pi t ) cos(2 \pi x)
     */
    
//     Term in the semidiscrete equation
    double semiDisc = 0;
    
//    These are set with boundary conditions
    int nextSpaceInd = 0;
    int prevSpaceInd = 0;
    
    int numTimeSteps = int(this->data.size1());
    int numSpaceSteps = int(this->data.size2());

    for (int timeStep = 1; timeStep < numTimeSteps; timeStep++) {
        
        for (int xStep = 0; xStep < numSpaceSteps; xStep++) {
            if (xStep == 0) {
                prevSpaceInd = int(numSpaceSteps) - 1; // for the k-1 space contribution at the beginning of the row
                semiDisc = data(timeStep - 1, xStep + 1) - 2 * data(timeStep -1, xStep) + data(timeStep - 1, prevSpaceInd);
            } else if (xStep == numSpaceSteps - 1) {
                nextSpaceInd = 0; // for the k+1 space contribution at the end of the row
                semiDisc = data(timeStep - 1, nextSpaceInd) - 2 * data(timeStep -1, xStep) + data(timeStep - 1, xStep - 1);
            } else {
                semiDisc = data(timeStep - 1, xStep + 1) - 2 * data(timeStep -1, xStep) + data(timeStep - 1, xStep - 1);
            }
            data(timeStep, xStep) = data(timeStep - 1, xStep) + dt / (std::pow(dx, 2)) * semiDisc;
        }
        
    }
}

void Models::findRMSError() {
    double error = 0;
    ublas::vector<double> exact(21);
    ublas::vector<double> est(21);
    int finalTStep = int(this->data.size1()) - 1;
    for (int col = 0; col < this->data.size2(); col++) {
        error += std::pow(exactSoln(finalTStep, col) - data(finalTStep,col), 2);
        exact(col) = exactSoln(finalTStep, col);
        est(col) = data(finalTStep,col);
    }
    
    this->rms = sqrt(error / this->data.size2());
}

inline double Models::exactSoln(int timeStep, int spaceStep) {
    // (spacestep * dx) = current value of x
    // (timestep * dt)  = current value of t
    return 1.0 + exp(-4 * std::pow(M_PI, 2) * this->duration) * cos(2 * M_PI * (spaceStep * dx));
    
}

