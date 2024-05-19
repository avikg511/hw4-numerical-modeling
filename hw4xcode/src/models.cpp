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

void Models::runForwardEulerMatrix() {
    /*
     Forward Euler Code with Matrices:
     $\frac{\partial C}{ \partial t} = \frac{ \partial^2 C }{\partial x^2} \\
                        = \frac{1}{\Delta x^2} (C^n_{k+1} -2C^n_{k} + C^n_{k-1})$ \\
     so $C^{n+1}_k = C^n_{k} + \frac{ \Delta t }{\Delta x^2} (C^n_{k+1} -2C^n_{k} + C^n_{k-1})$
     
     Each row of the matrix is size 20, and we're multiplying by a size 20 vector, so our matrix
        is ${\bf R}^{20 x 20}$. We want to implement boundary conditions with the k-1 and k+1 issues at
        k = 0 and k = 19 (the last element).
     */
    
    int numTimeSteps = int(this->data.size1());
    int numSpaceSteps = int(this->data.size2());
    
    ublas::matrix<double> FEMatrix = ublas::identity_matrix<double>(numSpaceSteps);
    
//    Scale the diagonal by -2 for the -2 * C_{k} term
    FEMatrix *= -2;
    
//    Boundary Conditions
    FEMatrix(0, FEMatrix.size2() - 1) = 1;
    FEMatrix(FEMatrix.size1() - 1, 0) = 1;
    
//    Add in the other diagonals
    for (int row = 0; row < FEMatrix.size1(); row++) {
        if (row == 0) {
//             For the first row, we want the C_{k+1} term to be +1
            FEMatrix(row, row + 1) = 1;
        } else if (row == (FEMatrix.size1() - 1)) {
//             For the last row, we want the C_{k-1} term to be +1
            FEMatrix(row, row - 1) = 1;
        } else {
//             In general, we want the C_{k-1} and C_{k+1} to be +1
            FEMatrix(row, row - 1) = 1;
            FEMatrix(row, row + 1) = 1;
        }
    }
    
//    Initializing matrix rows with arbitrary values
    ublas::matrix_row<ublas::matrix<double> > nMinus1Row = ublas::matrix_row<ublas::matrix<double>>(data, 0);
    
    for (int timeStep = 1; timeStep < numTimeSteps; timeStep++) {
        nMinus1Row = ublas::matrix_row<ublas::matrix<double>>(data, timeStep - 1);
        auto mult = ublas::prod(FEMatrix, nMinus1Row) * dt / (std::pow(dx, 2));
        ublas::row(data, timeStep) = nMinus1Row + mult;
    }
}

void Models::runImplicitEulerMatrix() {
    /*
     Implicit Euler Math:
        \vec{C^{n+1}} and \vec{C^{n}} refer to all C values in a vector at the respective time steps (n+1, n)
     
        \vec{C^{n+1}} = \vec{C^{n}} + \frac{\Delta t}{\Delta x^2} * M \vec{C^{n+1}}
        \implies I * \vec{C^{n+1}} - M \vec{C^{n+1}} = \vec{C^{n}}
        \implies (I - M) \vec{C^{n+1}} = \vec{C^{n}}
        \impiies \vec{C^{n+1}} = inv(I - M) \vec{C^{n}},
        where M is the matrix of 1, -2, 1 down the tridiagonal and 1's in the top right and bottom left corners.
     */
    
    int numTimeSteps = int(this->data.size1());
    int numSpaceSteps = int(this->data.size2());
    
    ublas::matrix<double> IEMatrix = ublas::identity_matrix<double>(numSpaceSteps);
    ublas::identity_matrix<double> eye(numSpaceSteps);
    
//    Scale the diagonal by -2 for the -2 * C_{k} term
    IEMatrix *= -2;
    
//    Boundary Conditions
    IEMatrix(0, IEMatrix.size2() - 1) = 1;
    IEMatrix(IEMatrix.size1() - 1, 0) = 1;
    
//    Add in the other diagonals
    for (int row = 0; row < IEMatrix.size1(); row++) {
        if (row == 0) {
//             For the first row, we want the C_{k+1} term to be +1
            IEMatrix(row, row + 1) = 1;
        } else if (row == (IEMatrix.size1() - 1)) {
//             For the last row, we want the C_{k-1} term to be +1
            IEMatrix(row, row - 1) = 1;
        } else {
//             In general, we want the C_{k-1} and C_{k+1} to be +1
            IEMatrix(row, row - 1) = 1;
            IEMatrix(row, row + 1) = 1;
        }
    }
    
//    Initializing matrix rows with arbitrary values
    ublas::matrix_row<ublas::matrix<double> > nMinus1Row = ublas::matrix_row<ublas::matrix<double>>(data, 0);
    
//    Inverting the Matrix
    ublas::matrix<double> IEOperator(numSpaceSteps, numSpaceSteps);
    if (!Models::InvertMatrix(eye - IEMatrix * this->dt / (std::pow(this->dx, 2)), IEOperator)) {
        return ;
    }
    
    for (int timeStep = 1; timeStep < numTimeSteps; timeStep++) {
        nMinus1Row = ublas::matrix_row<ublas::matrix<double>>(data, timeStep - 1);
        auto prod = ublas::prod(IEOperator, nMinus1Row);
        ublas::row(data, timeStep) = prod;
    }
}

int calcIndexMatsuno(int curSpaceStep, int dist, int numSpaceSteps) {
    if (curSpaceStep + dist < 0) {
        return numSpaceSteps - abs(dist) - 1;
    } else if (curSpaceStep + dist > (numSpaceSteps - 1)) {
        return curSpaceStep + dist - numSpaceSteps;
    } else {
        return curSpaceStep + dist;
    }
}

void Models::runMatsuno() {
    /*
     Matsuno Time Stepping Formulation:
     C_k^{n+1} = C_k + \sigma (C_{k+1} - 2C_k + C_{k-1}) + \sigma^2 (C_{k+2} -4C_{k+1} +6C_{k} - 4C_{k-1} +C_{k+2})
        where the RHS's terms are all in the nth timestep.
     
     semiDiscTerm1 = (C_{k+1} - 2C_k + C_{k-1})
     semiDiscTerm2 = (C_{k+2} -4C_{k+1} +6C_{k} - 4C_{k-1} +C_{k+2})
     
     \sigma = \frac{ \Delta t }{ \Delta x^2 }
     */
    
//     Term in the semidiscrete equation
    double semiDiscTerm1 = 0;
    double semiDiscTerm2 = 0;
    double sigma = dt / (std::pow(dx, 2));
    
    int numTimeSteps = int(this->data.size1());
    int numSpaceSteps = int(this->data.size2());
    
//    All space steps from time n-1 to calculate the Kth space step at the nth time step.
    double CKPlus2 = 0;
    double CKPlus1 = 0;
    double CK = 0;
    double CKMinus1 = 0;
    double CKMinus2 = 0;

    for (int timeStep = 1; timeStep < numTimeSteps; timeStep++) {
        for (int xStep = 0; xStep < numSpaceSteps; xStep++) {
//            calcIndexMatsuno provides the correct value for the matrix to adjust for the wrapping boundary conditions.
            CKPlus2 = data(timeStep - 1, calcIndexMatsuno(xStep, 0 + 2, numSpaceSteps));
            CKPlus1 = data(timeStep - 1, calcIndexMatsuno(xStep, 0 + 1, numSpaceSteps));
            CK = data(timeStep - 1, calcIndexMatsuno(xStep, 0, numSpaceSteps));
            CKMinus1 = data(timeStep - 1, calcIndexMatsuno(xStep, 0 - 1, numSpaceSteps));
            CKMinus2 = data(timeStep - 1, calcIndexMatsuno(xStep, 0 - 2, numSpaceSteps));
            
            semiDiscTerm1 = CKPlus1 - 2 * CK + CKMinus1;
            semiDiscTerm2 = CKPlus2 - 4 * CKPlus1 + 6 * CK - 4 * CKMinus1 + CKMinus2;
            
            data(timeStep, xStep) = data(timeStep - 1, xStep) + sigma * semiDiscTerm1 + std::pow(sigma, 2) * semiDiscTerm2;
        }
        
    }
}

void Models::findRMSError() {
    double error = 0;
    ublas::vector<double> exact(20);
    ublas::vector<double> est(20);
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

// Copied from: https://gist.github.com/lilac/2464434 and simplified the Templating because the rest of the code
// expects doubles.
bool Models::InvertMatrix(const ublas::matrix<double>& input, ublas::matrix<double>& inverse) {
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;
    // create a working copy of the input
    matrix<double> A(input);
    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());

    // perform LU-factorization
    int res = int(lu_factorize(A,pm));
    
    if (res != 0 ) {
        return false;
    }

    // create identity matrix of "inverse"
    inverse.assign(ublas::identity_matrix<double>(A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
}
