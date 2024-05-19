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


    // static double dt;  //comMod.dt
    // static double saveIncr;  //Number of steps between saving results, comMod.saveIncr
    // static double nTS;   //Number of timesteps, comMod.nTS




    //静态变量，每个出口相同
    static int ASCII;

    //非静态变量，具有出口特异性
    double flowEachTime = 0.0;
    double preFrom1DEachTime = 0.0;
    int q = 1;  // Global Solution Loop
    cvOneDOptions opts;
    std::string outputFileName;
    double currentTime;
    bool wasSet = false;

    // Pointer to the Model
    cvOneDModel *model;

    // Solve the blood flow problem
    void prepro(void);

	  void QuerryModelInformation(void);
        
    //the main solve part
    void GenerateSolution(void);
    void Nonlinear_iter(int step);

    // Result Output
    void postprocess_Text();
    void postprocess_VTK();
    void postprocess_VTK_XML3D_ONEFILE();
    void postprocess_VTK_XML3D_MULTIPLEFILES();

    //一些没用过的函数
    // Find Segment index given the ID
    int getSegmentIndex(int segID);
    // Get the solution;
    double GetSolution(int i, int j){return TotalSolution[i][j];}//IV 082103


 private:

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

    // Query Model, Allocate Memory, Set Initial Conditions
    static void CreateGlobalArrays(void);

    //initialize the solution, flow rate and area
    static void CalcInitProps(long subdomainID);
    
    //create MthSegmentModel and MthBranchModel if exists. Also specify inflow profile
    static void DefineMthModels(void);
    static void AddOneModel(cvOneDMthModelBase* model);


};

#endif
