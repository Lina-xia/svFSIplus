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

//
//  LinearSolver.cxx - Source for a Linear Skyline Matrix Solver
//  ~~~~~~~~~~~~~~~~
//
//  This class provides functionality for solving matrix systems
//  presented in the skyline format, and some special manipulations.
//


# include <cassert>

# include "cvOneDSkylineLinearSolver.h"
# include "cvOneDSkylineMatrix.h"
# include "cvOneDGlobal.h"

#define dmax(a,b) ((a<b) ? (b) : (a))


cvOneDSkylineLinearSolver::cvOneDSkylineLinearSolver(){

}

cvOneDSkylineLinearSolver::~cvOneDSkylineLinearSolver(){

}

void cvOneDSkylineLinearSolver::SetLHS(cvOneDFEAMatrix* matrix){
  lhsMatrix = matrix;
}

void cvOneDSkylineLinearSolver::SetRHS(cvOneDFEAVector *vector){
  rhsVector = vector;
}

void cvOneDSkylineLinearSolver::Solve(cvOneDFEAVector& sol){

  assert( rhsVector->GetDimension() == lhsMatrix->GetDimension());
  // cout<< *lhsMatrix;
  long numberOfEquations = lhsMatrix->GetDimension();
  long* position = ((cvOneDSkylineMatrix*)lhsMatrix)->GetPosition();
  double* KU = ((cvOneDSkylineMatrix*)lhsMatrix)->GetUpperDiagonalEntries();
  double* KL = ((cvOneDSkylineMatrix*)lhsMatrix)->GetLowerDiagonalEntries();
  double* KD = ((cvOneDSkylineMatrix*)lhsMatrix)->GetDiagonalEntries();
  double* F = rhsVector->GetEntries();
  double* solution = sol.GetEntries();
 //count size and nnz


  SolNonSymSysSkyLine( KU, KL, KD, F, position, solution, numberOfEquations, 0, EPSILON);
  SolNonSymSysSkyLine( KU, KL, KD, F, position, solution, numberOfEquations, 1, EPSILON);
}

cvOneDFEAMatrix* cvOneDSkylineLinearSolver::GetLHS(){
  return lhsMatrix;
}

cvOneDFEAVector* cvOneDSkylineLinearSolver::GetRHS(){
  return rhsVector;
}

void cvOneDSkylineLinearSolver::SetSolution(long equation, double value){

  // nent is the same for the corresponding column
  long nent = ((cvOneDSkylineMatrix*)lhsMatrix)->GetNumberOfEntriesIn(equation);

  double* columnValues = new double[nent];
  long* entries = new long[nent];
  assert( columnValues != 0 && entries != 0);

  ((cvOneDSkylineMatrix*)lhsMatrix)->GetColumnEntries( equation, entries, columnValues);

  lhsMatrix->ClearRow( equation);
  lhsMatrix->ClearColumn( equation);

  lhsMatrix->SetValue( equation, equation, 1.0);
  (*rhsVector)[equation] = value;

  // subtract values from the right hand side
  long i;
  //used?
  for( i = 0; i < nent; i++){
    (*rhsVector)[entries[i]] -= value * columnValues[i];
  }

  delete [] columnValues;
  delete [] entries;
}

// The values of one dense matrix are to be changed
//                                k_m-1,m^11 | k_m-1,m^12(kr0)    row 1
//                                k_m-1,m^21 | k_m-1,m^22(kr1)    row 2
//  k_m,m-1^11     k_m,m-1^12     k_m,m^11   | k_m,m^12(kr2)      row 3
//  ------------------------------------------
//  k_m,m-1^11(k0) k_m,m-1^12(k1) k_m,m^21(k2) k_m,m^22(k3)      row 4
//    col 1       col 2       col 3        col 4
//
//  At first, multiply row 4 by k_m, then add it to row 3;
//  secondly multiply column 4 by k_m, then add it to column 3;
//  finally the upper 3x3 matrix is kept.

void cvOneDSkylineLinearSolver::Minus1dof(long rbEqnNo, double k_m){
  double k[3]; //the array of k1...k3 shown above
  double kr[3]; //the array of kr1...kr3 shown above
  int i;
  double* KL = ((cvOneDSkylineMatrix*)lhsMatrix)->GetLowerDiagonalEntries();
  double* KU = ((cvOneDSkylineMatrix*)lhsMatrix)->GetUpperDiagonalEntries();
  double* KD = ((cvOneDSkylineMatrix*)lhsMatrix)->GetDiagonalEntries();
  for(i=3;i>0;i--){
    long pos = ((cvOneDSkylineMatrix*)lhsMatrix)->GetPosition(rbEqnNo, rbEqnNo-i);
    k[3-i] = KL[pos];
    kr[3-i] = KU[pos];
    KL[pos] = 0;
    KU[pos] = 0;
  }
  kr[2] += KD[rbEqnNo]*k_m;
  KD[rbEqnNo] = 1;
  for(i = 3; i > 0; i--){
    lhsMatrix->AddValue(rbEqnNo-1, rbEqnNo-i, k[3-i]*k_m);
    lhsMatrix->AddValue(rbEqnNo-i, rbEqnNo-1, kr[3-i]*k_m);
  }
  (*rhsVector)[rbEqnNo-1] += (*rhsVector)[rbEqnNo]*k_m;
  (*rhsVector)[rbEqnNo] = 0;
}

void cvOneDSkylineLinearSolver::DirectAppResistanceBC( long rbEqnNo, double resistance, double dpds, double rhs){
  int i;
  double* KL = ((cvOneDSkylineMatrix*)lhsMatrix)->GetLowerDiagonalEntries();
  double* KD = ((cvOneDSkylineMatrix*)lhsMatrix)->GetDiagonalEntries();
  long pos;

  for(i=3;i>1;i--){
    pos = ((cvOneDSkylineMatrix*)lhsMatrix)->GetPosition(rbEqnNo, rbEqnNo-i);
    KL[pos] = 0;
  }
  pos = ((cvOneDSkylineMatrix*)lhsMatrix)->GetPosition(rbEqnNo, rbEqnNo-1);
  KL[pos] = dpds;
  KD[rbEqnNo] = - resistance;
  (*rhsVector)[rbEqnNo] = rhs;
}

//assumes 2 nodes/element and 2degrees of freedom/node
void cvOneDSkylineLinearSolver::AddFlux(long rbEqnNo, double* OutletLHS11, double* OutletRHS1){

  lhsMatrix->AddValue(rbEqnNo-1, rbEqnNo-1, *OutletLHS11);
  lhsMatrix->AddValue(rbEqnNo-1, rbEqnNo, *(OutletLHS11+1));
  lhsMatrix->AddValue(rbEqnNo, rbEqnNo-1, *(OutletLHS11+2));
  lhsMatrix->AddValue(rbEqnNo, rbEqnNo, *(OutletLHS11+3));

  //cout<<rbEqnNo<<" "<<rbEqnNo-1<<endl;

  (*rhsVector)[rbEqnNo-1] += *OutletRHS1;
  (*rhsVector)[rbEqnNo] += *(OutletRHS1+1);
}

int cvOneDSkylineLinearSolver::SolNonSymSysSkyLine( double *KS, double *KI, double *KD, double *F,
                                             long *maxa, double *u, long neq, int solve,
                                             double eps){
  long it, firstp, fpvec, firstpscalp, len, i, posd1, posd2, pos1, pos2, pos, posd;
  double *vecLI, *vecCI, *vecLS, *vecCS, *vecLD, *vecCD;
  double sumI, sumS, sumD;

  //solves system. K has already been decomposed in LU form
  if(solve){
    solvLT( KI, F, maxa, neq); // solves L u'= F -> F := u'
    solvUT( KS, KD, u, F, maxa, neq);   // solves U u = u' = F
    return(1);
  }
  // Decomposes K in LU form.
  for( it = 1; it < neq; it++){

    firstp = it - (maxa[it+1] - maxa[it]);
    posd1 = maxa[it + 1];

    for( i = firstp; i < it; i++){

      fpvec = i - (maxa[i+1] - maxa[i]);
      firstpscalp = dmax( fpvec, firstp);
      len = (i - 1) - firstpscalp + 1;

      posd2 = maxa[ i + 1];
      pos1  = it - firstpscalp;
      pos2  =  i - firstpscalp;

      vecLI = &KI[posd1 - pos1];
      vecCI = &KS[posd2 - pos2];
      vecLS = &KI[posd2 - pos2];
      vecCS = &KS[posd1 - pos1];

      sumI = scalv( vecLI, vecCI, len);
      sumS = scalv( vecLS, vecCS, len);

      pos = posd1 - (it - i);
      KI[pos] = (KI[pos] - sumI) / KD[i];
      KS[pos] -= sumS;
    }

    len = (it - 1) - firstp + 1;
    posd = posd1 - (it - firstp);

    vecLD = &KI[posd];
    vecCD = &KS[posd];

    sumD = scalv( vecLD, vecCD, len);
    KD[it] -= sumD;

    //cout<<KD[it]<<endl;
    if (fabs(KD[it]) < eps){
      return(0);
    }
  }
  return(1);
}

void cvOneDSkylineLinearSolver::solvLT(double *KI, double *F, long *maxa, long neq){

  long i, firstp, len;
  double *vecA, *vecB;
  double sum;

  for( i = 0; i < neq; i++){
    firstp = i - (maxa[i+1] - maxa[i]);
    vecA = &KI[maxa[i]];
    vecB =  &F[firstp];
    len = (i - 1) - firstp + 1;
    sum = scalv( vecA, vecB, len);
    F[i] -= sum;
  }
}

void cvOneDSkylineLinearSolver::solvUT(double *KS, double *KD, double *u, double *F,
                                long *maxa, long neq){
  long i, j, it, firstp;

  u[neq - 1] = F[neq - 1] / KD[neq - 1];

  for( i = neq - 2; i > -1; i--){
    it = i + 1;
    firstp = it - (maxa[it + 1] - maxa[it]);
    for( j = firstp; j < it; j++){
      F[j] -= u[it] * KS[maxa[it + 1] - it + j];
    }
    u[i] = F[i] / KD[i];
  }
}

double cvOneDSkylineLinearSolver::scalv(double *p1, double *p2, long len){
  double sum = 0.0;
  while(len--){
    sum += (*p1++) * (*p2++);
  }
  return(sum);
}
