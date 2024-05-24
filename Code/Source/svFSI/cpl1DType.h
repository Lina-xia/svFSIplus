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
  public:

    // static double dt;  //comMod.dt
    // static double saveIncr;  //Number of steps between saving results, comMod.saveIncr
    // static double nTS;   //Number of timesteps, comMod.nTS

    //静态变量表示所有的出口共享这些数据
    static int ASCII;
    static std::string OutputFile;
    static bool path;

    //非静态变量，具有出口特异性
    std::string outletName;
    std::string inputFileName;

    double flowEachTime = 0.0;
    double preFrom1DEachTime = 0.0;
    int    q = 1;  // Global Solution Loop
    double currentTime;
    bool   wasSet = false;

    cvOneDOptions opts;
    // Pointer to the Model
    cvOneDModel *model;

    // Solve the blood flow problem
    void prepro(void);
	  void QuerryModelInformation(void);
    void GenerateSolution(void);
    void Nonlinear_iter(int step);

    // Result Output
    void postprocess_Text(std::string& path);
    void postprocess_VTK_XML3D_ONEFILE(std::string& path);
    void postprocess_VTK_XML3D_MULTIPLEFILES(std::string& path);

    //一些没用过的函数
    // Find Segment index given the ID
    int getSegmentIndex(int segID);
    // Get the solution;
    double GetSolution(int i, int j){return TotalSolution[i][j];}//IV 082103


 private:

    vector<cvOneDSubdomain*> subdomainList;
    vector<cvOneDFEAJoint*> jointList;
    vector<int> outletList;
    cvOneDFEAVector *currentSolution;
    cvOneDFEAVector *previousSolution;
    cvOneDFEAVector *increment;
    cvOneDFEAVector *rhs;
    cvOneDFEAVector *relLength; //for pressure calculation
    cvOneDMatrix<double> TotalSolution;
    cvOneDFEAMatrix *lhs;  // Generic matrix, can be skyline or sparse
    vector<cvOneDMthModelBase*> mathModels;

    // Query Model, Allocate Memory, Set Initial Conditions
    void CreateGlobalArrays(void);

    //initialize the solution, flow rate and area
    void CalcInitProps(long subdomainID);
    
    //create MthSegmentModel and MthBranchModel if exists. Also specify inflow profile
    void DefineMthModels(void);
    void AddOneModel(cvOneDMthModelBase* model);


};

#endif
