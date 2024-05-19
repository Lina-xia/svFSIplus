#ifndef CPL1DTYPE_H
#define CPL1DTYPE_H

# include <vector>
# include <iostream>
# include <ostream>
# include <cstdlib>
# include <cstdio>
# include <math.h>

# include "cvOneDTypes.h"
# include "cvOneDModel.h"
# include "cvOneDSubdomain.h"
# include "cvOneDMthModelBase.h"
# include "cvOneDFEAJoint.h"
# include "cvOneDOptions.h"


class cpl1DType
{
  public://静态变量表示所有的出口共享这些数据

    double flowEachTime = 0.0;
    double preFrom1DEachTime = 0.0;
    // Global Solution Loop
    int q = 1;

    // static double dt;  //comMod.dt
    // static double saveIncr;  //Number of steps between saving results, comMod.saveIncr
    // static double nTS;   //Number of timesteps, comMod.nTS
    
    // couple1D
    cvOneDOptions opts;
    std::string outputFileName;

    static double currentTime;
    static int ASCII;

    // Set the Model Pointer
    static void SetModelPtr(cvOneDModel *mdl);
    static cvOneDModel* GetModelPtr(){return model;}

    // Solve the blood flow problem
    static void Solve(void);

    // Get the solution;
    static double GetSolution(int i, int j){return TotalSolution[i][j];}//IV 082103

    // Cleanup
    static void Cleanup(void);

    static void DefineInletFlow(double* time, double* flrt, int num);

	  static void QuerryModelInformation(void);


    // Result Output
    static void postprocess_Text();
    static void postprocess_VTK();
    static void postprocess_VTK_XML3D_ONEFILE();
    static void postprocess_VTK_XML3D_MULTIPLEFILES();

    // Find Segment index given the ID
    static int getSegmentIndex(int segID);
    //the main solve part
    static void GenerateSolution(void);
    void Nonlinear_iter(int step);

 private:

    // Query Model, Allocate Memory, Set Initial Conditions
    static void CreateGlobalArrays(void);

    //initialize the solution, flow rate and area
    static void CalcInitProps(long subdomainID);
    
    //create MthSegmentModel and MthBranchModel if exists. Also specify inflow profile
    static void DefineMthModels(void);
    static void AddOneModel(cvOneDMthModelBase* model);

    static bool wasSet;

    // Pointer to the Model
    static cvOneDModel *model;

    static vector<cvOneDSubdomain*> subdomainList;
    static vector<cvOneDFEAJoint*> jointList;
    static vector<int> outletList;
    static cvOneDFEAVector *currentSolution;
    static cvOneDFEAVector *previousSolution;
    static cvOneDFEAVector *increment;
    static cvOneDFEAVector *rhs;
    static cvOneDFEAVector *relLength; //for pressure calculation
    static cvOneDMatrix<double> TotalSolution;

    // Generic matrix, can be skyline or sparse
    static cvOneDFEAMatrix *lhs;

    static vector<cvOneDMthModelBase*> mathModels;

    static long numFlowPts;
    static double *flowRate;
    static double *flowTime;
	  static double Period;

};

#endif
