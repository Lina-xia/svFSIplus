/* Copyright (c) Stanford University, The Regents of the University of
 *               California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef CVONEDLINEARSOLVER_H
#define CVONEDLINEARSOLVER_H

//
//  cvOneDLinearSolver.h - Header for a Linear Skyline Matrix Solver
//  
//  This class provides functionality for solving matrix systems
//  presented in the skyline format, and some special manipulations. 
//  

# include <cmath>

# include "cvOneDFEAMatrix.h"
# include "cvOneDFEAVector.h"

class cvOneDLinearSolver{

  public:

    cvOneDFEAMatrix* lhsMatrix;
    cvOneDFEAVector* rhsVector;

    cvOneDLinearSolver();
    virtual ~cvOneDLinearSolver();

    virtual void SetLHS( cvOneDFEAMatrix* matrix) = 0;  //有替代
    virtual void SetRHS( cvOneDFEAVector* vector) = 0;  //有替代
  
    // matrix is overwritten with its LU decomposition
    // solution gets overwritten with the solution of 
    // the linear system of equations
    virtual void Solve(cvOneDFEAVector& solution) = 0;
  
    virtual cvOneDFEAMatrix* GetLHS() = 0;
    virtual cvOneDFEAVector* GetRHS() = 0;
    virtual void SetSolution(long equation, double value) = 0;
    // when one more constraint (dQ = k_m*dS, resistance boundary
    // condition) is added, the basic dense matrix (4x4) is decreased
    // to (3x3)
    virtual void Minus1dof(long rightBottomEquationNumber, double k_m) = 0;
    // direct application of the resistance constraint without reduction, which means 
    // Newton Raphson scheme is applied on equation Q = PR. Therefore the dense matrix 
    // is still 4x4. 
    virtual void DirectAppResistanceBC(long rbEqnNo, double resistance, double dpds, double rhs) = 0;
	  virtual void AddFlux(long rbEqnNo, double* OutletLHS11, double* OutletRHS1) = 0;
	
};

#endif // CVONEDLINEARSOLVER_H
