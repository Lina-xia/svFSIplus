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

using namespace std;

class cpl1DType
{
  public:

    //静态变量表示所有的出口共享这些数据
    static bool     IfCout;
    static double   dt;  //comMod.dt
    static int      saveIncr;  //comMod.saveIncr
    static int      maxStep;   //comMod.nTS
    static string   OutputFile;
    // SOLVER OPTIONS
    static int    quadPoints;
    static double convergenceTolerance;
    static int    useIV;
    static int    useStab;
    static bool   solverOptionDefined;
    // MATERIAL
    static vector<string> materialName;
    static vector<string> materialType;
    static vector<double> materialDensity;
    static vector<double> materialViscosity;
    static vector<double> materialPRef;
    static vector<double> materialExponent;
    static vector<double> materialParam1;
    static vector<double> materialParam2;
    static vector<double> materialParam3;
    //OUTPUT
    static int    outputType;    //0表示TEXT，1表示VTP
    static int    vtkOutputType;
    static int    ASCII;         //VTK输出格式


    //非静态变量，具有出口特异性
    cvOneDOptions opts;
    cvOneDModel   *model;  // Pointer to the Model
    string   outletName;
    string   inputFile;
    double        flowEachTime = 0.0;
    double        preFrom1DEachTime = 0.0;
    int           q = 1;  // Global Solution Loop
    bool          wasSet = false;

    // Solve the blood flow problem
    void readSharedVar();
    void prepro(void);
	  void QuerryModelInformation(void);
    void GenerateSolution(void);
    void Nonlinear_iter(int step);
    int  getDataTableIDFromStringKey(string key);
    void createModel();
    void readModelFile(cvStringVec includedFiles);
    void readModel();

    // Result Output
    void postprocess_Text(std::string& path);
    void postprocess_VTK_XML3D_ONEFILE(std::string& path);
    void postprocess_VTK_XML3D_MULTIPLEFILES(std::string& path);
    
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
