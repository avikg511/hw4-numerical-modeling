//
//  main.cpp
//  hw4xcode
//
//  Created by Avik Ghosh on 5/11/24.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include "models.hpp"
#include <fstream>

int main(int argc, const char * argv[]) {
    /*
     As shown in the writeup, $\sigma = \frac{\Delta t}{\Delta x^2}$
     
     We want to show that $\sigma = \frac{1}{6}$ gives minimum error, and given that we have $\Delta x = 0.05,
        we can figure out what $\Delta t$ should be now. Multiplying both sides by $0.05^2$, we see that
        the optimal $\Delta t$, $dt^*$ = 0.00004464285714 $\approx$ 4.464E-5. That means we still need a pretty fine
        resolution.
     
     The way we'll pick the range is basically $\sigma \in [1E-7, 1)$ with 5E-7 resolution. Then scale it
        by 0.05^2 to designate them all as $\Delta t$ values. In the end, we'll plot $\Delta t$ on the x axis.
     
     */
    long double duration = 0.3;
    long double xStepSize = 0.05;
    
    // Sigma values here are the variables we modify to test numerical stability/rms error
    ublas::vector<long double> sigma (1800);
    boost::algorithm::iota(sigma.begin(), sigma.end(), 1);

    // Currently stepped through from 1 to 1800, and then scaled by time/space resolution
    sigma *= 5;
    sigma *= 1.0 / std::pow(10, 7) / std::pow(xStepSize, 2);

    // Time Values here for timesteps, it's just the sigma value with the space resolution multiplied
    auto tVals = sigma * std::pow(xStepSize, 2);
    
    // Data outputs to this file
    std::ofstream myFile("outputData.csv");
    
    for (auto it = tVals.begin(); it != tVals.end(); it++) {
//         It's not my favorite piece of code, but here I just reinitialize the class every timestep. Otherwise, \\
            I'd need to basically delete the matrices, reinitialize the sizing with the new time resolution, etc. etc.
        Models model(*it, duration, xStepSize);
        model.runImplicitEulerMatrix(); // Switch out with different methods as needed
        model.findRMSError();
        myFile << *it << "," << model.rms << std::endl;

    }
    return 0;
}
