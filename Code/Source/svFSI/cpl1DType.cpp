///////////////////////
// c p l 1 D T y p e //
///////////////////////

# include <time.h>
# include <iomanip> // 包含 std::setprecision 和 std::scientific 等操纵器

# include "cpl1DType.h"
# include "cvOneDGlobal.h"
# include "cvOneDString.h"
# include "cvOneDException.h"
# include "cvOneDMaterial.h"
# include "cvOneDMthSegmentModel.h"
# include "cvOneDMthBranchModel.h"
# include "cvOneDModelManager.h"

#ifndef WIN32
#define _USE_MATH_DEFINES
#endif

# ifdef USE_SKYLINE
  # include "cvOneDSkylineMatrix.h"
  # include "cvOneDSkylineLinearSolver.h"
# endif

# ifdef USE_SUPERLU
  # include "sparse/cvOneDSparseMatrix.h"
  # include "sparse/cvOneDSparseLinearSolver.h"
# endif

# ifdef USE_CSPARSE
  # include "sparse/cvOneDSparseMatrix.h"
  # include "sparse/cvOneDSparseLinearSolver.h"
# endif

# define baryeTommHg 0.0007500615613026439

using namespace std;

// Static Declarations...
bool      cpl1DType::CreateCout = true;
bool      cpl1DType::solverOptionDefined = false;
// 需要广播的数据
int       cpl1DType::quadPoints = 2;
double    cpl1DType::convergenceTolerance = 1.0e-8;
int       cpl1DType::useIV = 1;
int       cpl1DType::useStab = 1;
int       cpl1DType::outputType = 0; // Default Text Output; 0表示TXET;1表示VTP;2表示BOTH
string    cpl1DType::OutputFile = string("OneDSolver.out");

vector<string> cpl1DType::materialName = {};
vector<string> cpl1DType::materialType = {};
vector<double> cpl1DType::materialDensity;
vector<double> cpl1DType::materialViscosity;
vector<double> cpl1DType::materialPRef;
vector<double> cpl1DType::materialExponent;
vector<double> cpl1DType::materialParam1;
vector<double> cpl1DType::materialParam2;
vector<double> cpl1DType::materialParam3;


// ==================
// WRITE TEXT RESULTS
// ==================
void cpl1DType::postprocess_Text(string& path){
  int j;
  int fileIter;
  int elCount = 0;
  // cout << model->getModelName() << endl;
  // cout << model->getModelID() <<endl;
  for(fileIter = 0; fileIter < model -> getNumberOfSegments(); fileIter++){

    cvOneDSegment *curSeg = model -> getSegment(fileIter);
    cvOneDMaterial* curMat = subdomainList[fileIter]->GetMaterial();

    long numEls = curSeg -> getNumElements();
    double segLength = curSeg->getSegmentLength();

    long startOut = elCount;
    long finishOut = elCount + ((numEls+1)*2);

    char *tmp1 = curSeg -> getSegmentName();

    char tmp2[512];
    char tmp3[512];
    char tmp4[512];
    char tmp5[512];
    char tmp6[512]; // WSS

    strcpy(tmp2, path.c_str());
    strcpy(tmp3, path.c_str());
    strcpy(tmp4, path.c_str());
    strcpy(tmp5, path.c_str());
    strcpy(tmp6, path.c_str());

    char *btemp= model-> getModelName(); // add to write out binary files for java
    
    strcat(tmp2, btemp);
    strcat(tmp3, btemp);
    strcat(tmp4, btemp);
    strcat(tmp5, btemp);
    strcat(tmp6, btemp);

    strcat(tmp2, tmp1);
    strcat(tmp2, "_flow.csv");
    // cout << tmp2 << endl;
    strcat(tmp3, tmp1);
    strcat(tmp3, "_area.csv");
    strcat(tmp4, tmp1);
    strcat(tmp4, "_pressure.csv");
    // strcat(tmp5, tmp1);
    // strcat(tmp5,"_Re.csv");
    // strcat(tmp6, tmp1);
    // strcat(tmp6,"_wss.csv");

    ofstream flow,area,pressure,reynolds,wss,vecolity_flux; // for ASCII

    if(CreateFile){
      flow.open(tmp2, ios::out);
      area.open(tmp3, ios::out);
      pressure.open(tmp4, ios::out);
      // reynolds.open(tmp5, ios::out);
      // wss.open(tmp6, ios::out);
      CreateFile = false;
    }else{
      // Text Files
      flow.open(tmp2, ios::app);
      area.open(tmp3, ios::app);
      pressure.open(tmp4, ios::app);
      // reynolds.open(tmp5, ios::app);
      // wss.open(tmp6, ios::app);
    }
    flow.precision(OUTPUT_PRECISION);
    area.precision(OUTPUT_PRECISION);
    pressure.precision(OUTPUT_PRECISION);
    // reynolds.precision(OUTPUT_PRECISION);
    // wss.precision(OUTPUT_PRECISION);

    double val;

    // Output the flow file
    for(j=startOut+1;j<finishOut;j+=2){
        int i = step; //最后一个时间步的
        flow << TotalSolution[i][j] << ",";
    }


    // Output the Area/Pressure/Reynolds/WSS file
    int ii;
    double Area, Re, Pre, r, flo, wssVal,vel_flux;

    for(ii=0,j=startOut;ii<numEls+1 && j<finishOut;ii++,j+=2){
      double z = (ii/double(numEls))*segLength;
      int i= step;
        //write area
        Area = (double)TotalSolution[i][j];
        area << Area << ",";

        // Write Re
        r = sqrt(Area/M_PI);
        flo = (double)TotalSolution[i][j+1];
        //usual Re=rho/mu*D*velocity=rho/mu*Q*sqrt(4/Pi/Area)
        Re = curMat->GetDensity()/curMat->GetDynamicViscosity()*flo/sqrt(Area)*sqrt(4.0/M_PI);
        // reynolds << Re << ",";

        // Write pressure - Initial (CGS) Units
        Pre = (double) curMat->GetPressure(TotalSolution[i][j],z);
        pressure << curMat->GetPressure(TotalSolution[i][j],z) << ",";

        // // write Wall Shear Stresses for Poiseuille flow
        // wssVal = (4.0*curMat->GetDynamicViscosity()*flo)/(M_PI*r*r*r);
        //  wss << wssVal << ",";
    }

    elCount += 2*(numEls+1);

    area << endl;
    pressure << endl;
    flow << endl;
    // reynolds << endl;
    // wss << endl;
    area.close();
    pressure.close();
    flow.close();
    // reynolds.close();
    // wss.close();
  }
}

// Get the index of a segment given its ID
int cpl1DType::getSegmentIndex(int segID){
  cvOneDSegment* curSeg;
  bool found = false;
  int count = 0;
  while((!found)&&(count<model->getNumberOfSegments())){
    curSeg = model->getSegment(count);
    found = (curSeg->getSegmentID() == segID);
    if(!found){
      count++;
    }
  }
  if(found){
    return count;
  }else{
    throw cvException("ERROR:Cannot Find Segment ID in cpl1DType::getSegmentIndex\n");
  }
}

// EVAL SEGMENT AXIS SYSTEM
void evalSegmentLocalAxis(double axis[][3]){
  // Get Reference Axis Y or Z
  double ref[3];
  double mod = 0.0;
  double yDist = sqrt((axis[0][0] - 0.0)*(axis[0][0] - 0.0) +
                      (axis[1][0] - 1.0)*(axis[1][0] - 1.0) +
                      (axis[2][0] - 0.0)*(axis[2][0] - 0.0));
  if(yDist < 1.0e-5){
    // Reference Axis Z
    ref[0] = 0.0;
    ref[1] = 0.0;
    ref[2] = 1.0;
  }else{
    // Reference Axis Y
    ref[0] = 0.0;
    ref[1] = 1.0;
    ref[2] = 0.0;
  }
  // Make the First Outer Product
  axis[0][1] =  axis[1][0]*ref[2] - ref[1]*axis[2][0];
  axis[1][1] = -axis[0][0]*ref[2] + ref[0]*axis[2][0];
  axis[2][1] =  axis[0][0]*ref[1] - ref[0]*axis[1][0];
  mod = sqrt(axis[0][1]*axis[0][1] + axis[1][1]*axis[1][1] + axis[2][1]*axis[2][1]);
  axis[0][1] /= mod;
  axis[1][1] /= mod;
  axis[2][1] /= mod;
  // Make the Second Outer Product
  axis[0][2] =  axis[1][0]*axis[2][1] - axis[2][0]*axis[1][1];
  axis[1][2] = -axis[0][0]*axis[2][1] + axis[2][0]*axis[0][1];
  axis[2][2] =  axis[0][0]*axis[1][1] - axis[1][0]*axis[0][1];
  mod = sqrt(axis[0][2]*axis[0][2] + axis[1][2]*axis[1][2] + axis[2][2]*axis[2][2]);
  axis[0][2] /= mod;
  axis[1][2] /= mod;
  axis[2][2] /= mod;
}

// =======================================================
// WRITE 3D XML VTK RESULTS - MULTIPLE FILES FOR ANIMATION
// =======================================================
void cpl1DType::postprocess_VTK(std::string& path,int& cTS, double& scF){
  // Set Constant Number of Subdivisions on the vessel circumference
  int circSubdiv = 20;

  int currSegID = 0;
  cvOneDSegment* currSeg = NULL;
  cvOneDJoint* currJoint = NULL;
  cvOneDMaterial* curMat = NULL;
  cvOneDNode* currNode = NULL;

  // DEFINE INCIDENCE
  std::vector<double> segInlets(model->getNumberOfSegments());
  std::vector<double> segOutlets(model->getNumberOfSegments());
  for(int loopSegment=0;loopSegment<model->getNumberOfSegments();loopSegment++){
    segInlets[loopSegment] = -1;
    segOutlets[loopSegment] = -1;
  }

  // FORM INCIDENCE AND STORE COORDS
  cvDoubleMat nodeList;
  cvDoubleVec temp;
  long* segNodes;
  int inletNodeID = 0;
  int outletNodeID = 0;
  if(model->getNumberOfNodes() > 0){
    for(int loopNode = 0; loopNode < model->getNumberOfNodes(); loopNode++){
      temp.clear();
      // Get joint coordinate
      currNode = model->getNode(loopNode);
      temp.push_back(currNode->x);
      temp.push_back(currNode->y);
      temp.push_back(currNode->z);
      nodeList.push_back(temp);
    }
    // Loop on the segments
    for(int loopSegment = 0; loopSegment < model->getNumberOfSegments(); loopSegment++){
      // Get current segment
      currSeg = model->getSegment(loopSegment);
      // Get End Nodes
      segNodes = currSeg->getInOutJoints();
      // Mark Inlets
      inletNodeID = segNodes[0];
      outletNodeID = segNodes[1];
      // Assign Inlets and Outlets
      segInlets[loopSegment] = inletNodeID;
      segOutlets[loopSegment] = outletNodeID;
    }
  }else{
    // There are not Joints and a single segment
    currSeg = model->getSegment(0);
    // First Node
    temp.clear();
    temp.push_back(0.0);
    temp.push_back(0.0);
    temp.push_back(0.0);
    nodeList.push_back(temp);
    // Second Node
    temp.clear();
    temp.push_back(currSeg->getSegmentLength());
    temp.push_back(0.0);
    temp.push_back(0.0);
    nodeList.push_back(temp);
    // Store the connectivity
    segInlets[0] = 0;
    segOutlets[0] = 1;
  }

  for(int loopSegment=0;loopSegment<model->getNumberOfSegments();loopSegment++){
    if(segInlets[loopSegment] == -1){
      // CHECK INLETS
      printf("ERROR: INLET FOR SEGMENT %d\n",loopSegment);
    }
    if(segOutlets[loopSegment] == -1){
      // CHECK OUTLETS
      printf("ERROR: OUTLET FOR SEGMENT %d\n",loopSegment);
    }
  }

  // LOOP OVER TIME
  int totSegmentPoints = 0;
  int totSegmentSolutions = 0;
  long segOffset = 0;
  int inletSegJoint = 0;
  int outletSegJoint = 0;
  double segVers[3][3];
  double mod = 0.0;
  double currCentre[3] = {0.0};
  double currIniArea = 0.0;
  double currIniRad = 0.0;
  double currTheta = 0.0;
  cvDoubleVec tmp;
  cvDoubleMat segNodeList;
  double lengthByNodes = 0.0;
  double lengthBySegment = 0.0;
  int startOut = 0;
  int finishOut = 0;
  double segLength = 0.0;
  double z = 0.0;
  double flo = 0.0;
  double area = 0.0;
  double radius = 0.0;
  double Re = 0.0;
  double wss = 0.0;
  double iniArea = 0.0;
  double newArea = 0.0;
  double radDisp = 0.0;
  int loopTime = step; //只输出最后一个虚拟时间步的数据
  // cout << "step = " << step << endl;
  cvStringVec fileList;
  string fileName;

    // Set and open VTK file for current time step
    fileName = model->getModelName();
    char timeString[512];
    char suffix[512];
    sprintf(timeString, "_%05d", cTS);
    fileName = path.c_str() + fileName + string(timeString) + ".vtp";
    FILE* vtkFile;
    vtkFile = fopen(fileName.c_str(),"w");
    // Add to a list of files
    fileList.push_back(fileName);

    // Write VTK XML Header
    fprintf(vtkFile,"<?xml version=\"1.0\"?>\n");
    fprintf(vtkFile,"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(vtkFile,"<PolyData>\n");

    // Init Segment Offset
    segOffset = 0;

    // LOOP OVER THE SEGMENTS
    for(int loopSegment=0;loopSegment<model->getNumberOfSegments();loopSegment++){

      // Clear the node list for the segment
      segNodeList.clear();

      // Get Current Segment
      currSeg = model->getSegment(loopSegment);

      // Compute the total number of points for this segment
      totSegmentSolutions = (currSeg->getNumElements()+1);
      totSegmentPoints = totSegmentSolutions * circSubdiv;

      // Set the range for the totalsoluton of this segment
      startOut = segOffset;
      finishOut = segOffset + 2*(totSegmentSolutions);

      // Get Material
      curMat = subdomainList[loopSegment]->GetMaterial();

      // Write Piece Header
      fprintf(vtkFile,"<Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"%ld\" NumberOfPolys=\"0\">\n",totSegmentPoints,currSeg->getNumElements());

      // Get inlet and outlet joints
      inletSegJoint = segInlets[loopSegment];
      outletSegJoint = segOutlets[loopSegment];

      // Compute Segment Versor
      segVers[0][0] = nodeList[outletSegJoint][0] - nodeList[inletSegJoint][0];
      segVers[1][0] = nodeList[outletSegJoint][1] - nodeList[inletSegJoint][1];
      segVers[2][0] = nodeList[outletSegJoint][2] - nodeList[inletSegJoint][2];
      mod = sqrt(segVers[0][0]*segVers[0][0] + segVers[1][0]*segVers[1][0] + segVers[2][0]*segVers[2][0]);
      segVers[0][0] /= mod;
      segVers[1][0] /= mod;
      segVers[2][0] /= mod;

      lengthByNodes = mod;
      lengthBySegment = currSeg->getSegmentLength();

      // Compute Segment Local axis system
      evalSegmentLocalAxis(segVers);

      // Loop on the number of elements
      for(int loopEl=0;loopEl<currSeg->getNumElements() + 1;loopEl++){
        // Compress/Elongate solution by length between nodes rather than defined segment length - to prioritize segment length
        // and maintain similar geometry to node definitions, MD 4/2/19
        currCentre[0] = nodeList[inletSegJoint][0] + loopEl*lengthByNodes/double(currSeg->getNumElements())*segVers[0][0];
        currCentre[1] = nodeList[inletSegJoint][1] + loopEl*lengthByNodes/double(currSeg->getNumElements())*segVers[1][0];
        currCentre[2] = nodeList[inletSegJoint][2] + loopEl*lengthByNodes/double(currSeg->getNumElements())*segVers[2][0];

        // Get initial radius at current location
        currIniArea = currSeg->getInitInletS() + (loopEl/double(currSeg->getNumElements()))*(currSeg->getInitOutletS() - currSeg->getInitInletS());
        currIniRad = sqrt(currIniArea/M_PI);

        // Loop on the subdivisions
        for(int loopSubdiv=0;loopSubdiv<circSubdiv;loopSubdiv++){
          currTheta = loopSubdiv*2*M_PI/double(circSubdiv);
          tmp.clear();
          tmp.push_back(currCentre[0] + currIniRad*segVers[0][1]*cos(currTheta) + currIniRad*segVers[0][2]*sin(currTheta));
          tmp.push_back(currCentre[1] + currIniRad*segVers[1][1]*cos(currTheta) + currIniRad*segVers[1][2]*sin(currTheta));
          tmp.push_back(currCentre[2] + currIniRad*segVers[2][1]*cos(currTheta) + currIniRad*segVers[2][2]*sin(currTheta));
          segNodeList.push_back(tmp);
        }
      }

      // List of Node Coordinates ready for export
      fprintf(vtkFile,"<Points>\n");
      fprintf(vtkFile,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
      for(int loopA=0;loopA<segNodeList.size();loopA++){
        fprintf(vtkFile,"%e %e %e\n",segNodeList[loopA][0]/scF,segNodeList[loopA][1]/scF,segNodeList[loopA][2]/scF);
      }
      fprintf(vtkFile,"</DataArray>\n");
      fprintf(vtkFile,"</Points>\n");

      // Write Strip Incidence and offset
      fprintf(vtkFile,"<Strips>\n");
      // Strip Connectivity
      fprintf(vtkFile,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
      for(int loopA=0;loopA<currSeg->getNumElements();loopA++){
        for(int loopB=0;loopB<circSubdiv;loopB++){
          fprintf(vtkFile,"%d %d ",loopA*circSubdiv+loopB,loopA*circSubdiv+loopB+circSubdiv);
        }
        fprintf(vtkFile,"%d %d ",loopA*circSubdiv+0,loopA*circSubdiv+0+circSubdiv);
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");
      // Strip Offset
      fprintf(vtkFile,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
      for(int loopA=0;loopA<currSeg->getNumElements();loopA++){
        fprintf(vtkFile,"%d ",(loopA+1)*(circSubdiv*2+2));
      }
      fprintf(vtkFile,"\n");
      fprintf(vtkFile,"</DataArray>\n");
      fprintf(vtkFile,"</Strips>\n");


      // PRINT OUTPUTS
      fprintf(vtkFile,"<PointData Scalars=\"ScalOutputs\" Vectors=\"VecOutputs\">\n");

      // PRINT FLOW RATES
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Flowrate\" NumberOfComponents=\"1\" format=\"ascii\">\n");
      for(int j=startOut+1;j<finishOut;j+=2){
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",(double)TotalSolution[loopTime][j]);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT AREA
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Area\" NumberOfComponents=\"1\" format=\"ascii\">\n");
      for(int j=startOut;j<finishOut;j+=2){
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",(double)TotalSolution[loopTime][j]);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT RADIAL DISPLACEMENTS AS VECTORS
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Disps\" NumberOfComponents=\"3\" format=\"ascii\">\n");
      for(int j=startOut;j<finishOut;j+=2){

        // Evaluate Initial Area at current location
        iniArea = currSeg->getInitInletS() + (((j-startOut)/2)/double(currSeg->getNumElements()))*(currSeg->getInitOutletS() - currSeg->getInitInletS());
        // Eval Current Area at current location
        newArea = TotalSolution[loopTime][j];
        // Evaluate Radial displacement
        radDisp = sqrt(newArea/M_PI) - sqrt(iniArea/M_PI);
        for(int k=0;k<circSubdiv;k++){
          // Print the three components for every point
          currTheta = k*2*M_PI/double(circSubdiv);
          tmp.clear();
          tmp.push_back(radDisp*segVers[0][1]*cos(currTheta) + radDisp*segVers[0][2]*sin(currTheta));
          tmp.push_back(radDisp*segVers[1][1]*cos(currTheta) + radDisp*segVers[1][2]*sin(currTheta));
          tmp.push_back(radDisp*segVers[2][1]*cos(currTheta) + radDisp*segVers[2][2]*sin(currTheta));
          // Write values
          fprintf(vtkFile,"%e %e %e ",tmp[0],tmp[1],tmp[2]);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT PRESSURE IN MMHG
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Pressure_mmHg\" NumberOfComponents=\"1\" format=\"ascii\">\n");
      segLength = currSeg->getSegmentLength();
      int section = 0;
      for(int j=startOut;j<finishOut;j+=2){
        z = (section/(double)currSeg->getNumElements())*segLength;
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",curMat->GetPressure(TotalSolution[loopTime][j],z)*baryeTommHg);
        }
        fprintf(vtkFile,"\n");
        section++;
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT REYNOLDS NUMBER
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Reynolds\" NumberOfComponents=\"1\" format=\"ascii\">\n");
      for(int j=startOut;j<finishOut;j+=2){
        // Get Flow
        flo = (double)TotalSolution[loopTime][j+1];
        area = (double)TotalSolution[loopTime][j];
        // Get Area
        Re = curMat->GetDensity()/curMat->GetDynamicViscosity()*flo/sqrt(area)*sqrt(4.0/M_PI);
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",Re);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT WSS
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"WSS\" NumberOfComponents=\"1\" format=\"ascii\">\n");
      for(int j=startOut;j<finishOut;j+=2){
        // Get Flow
        flo = (double)TotalSolution[loopTime][j+1];
        // Get Radius
        radius = sqrt((double)TotalSolution[loopTime][j]/M_PI);
        // Get WSS
        wss = 4.0*curMat->GetDynamicViscosity()*flo/(M_PI*radius*radius*radius);
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",wss);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // Close Pointdata
      fprintf(vtkFile,"</PointData>\n");
      // Close Piece
      fprintf(vtkFile,"</Piece>\n");

      // Increment Segment Offset
      segOffset += 2*(totSegmentSolutions);

    } // End Segment Loop
    // Close
    fprintf(vtkFile,"</PolyData>\n");
    fprintf(vtkFile,"</VTKFile>\n");
    fclose(vtkFile);

}

// ====================
// MAIN SOLUTION DRIVER
// ====================
void cpl1DType::prepro(void){
  char errStr[256];
  // First check to make sure we've set a model pointer
  // Prior to solution attempt.
  if(model == NULL){
    cvOneDError::setErrorNumber(ErrorTypeScope::BAD_VALUE);
    strcpy(errStr,"In BFSolver::Solve(...), No model pointer was set prior to solution attempt.  Bailing out to avoid a coredump!");
    cvOneDError::setErrorString(errStr);
    cvOneDError::CallErrorHandler();
    exit(0);
  }

  // Query the model for information, this is where
  // the subdomain and material information get passed.
  QuerryModelInformation();

  DefineMthModels();

  CreateGlobalArrays();

  long id = model -> getNumberOfSegments();

  int i;
  for (i=0; i<id; i++){
    CalcInitProps(i);
  }

}


void cpl1DType::DefineMthModels(){
  mathModels.clear();

  // cout << "Subdomain No. "<<subdomainList.size() << endl;
  // cout << "Joint No. "<< jointList.size() << endl;
  // cout << "Outlet No. "<< outletList.size() << endl;
  cvOneDMthSegmentModel* segM = new cvOneDMthSegmentModel(subdomainList, jointList, outletList, quadPoints);

  cvOneDMthBranchModel* branchM = new cvOneDMthBranchModel(subdomainList, jointList, outletList);
  AddOneModel(segM);
  AddOneModel(branchM);
}

void cpl1DType::AddOneModel(cvOneDMthModelBase* model){
  mathModels.push_back(model);
}

void cpl1DType::QuerryModelInformation(void)
{
    // place to create the subdomain and material.
    long is, ij;
    is = model -> getNumberOfSegments();
    ij = model -> getNumberOfJoints();

    jointList.resize(0);
    outletList.resize(0);
    subdomainList.resize(0);

    // cout << endl;
    // cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    // cout << "Number of Joints: " << ij << endl;
    // cout << "Number of Segments: " << is << endl;
    // cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    // cout << endl;

    long int i, j;

    // Loop on Joints
    for(i=0; i<ij; i++){
      cvOneDFEAJoint *feajoint = new cvOneDFEAJoint();
      cvOneDJoint *joint = model->getJoint(i);
      for(j=0;j < model -> getJoint(i)->InletSegments.size(); j++){
        feajoint->AddInletSubdomains(joint->InletSegments[j]);
      }
      for(j=0;j < joint->OutletSegments.size(); j++){
        feajoint->AddOutletSubdomains(joint->OutletSegments[j]);
      }
      jointList.push_back(feajoint);
    }

    // Loop on Segments
    long temp = 0;
    for (i=0; i<is; i++){
      cvOneDSegment* seg = model->getSegment(i);
      long nels = seg->getNumElements();
      double segLen = seg->getSegmentLength();
      int matID = seg->getMaterialID();
      MeshType mType = seg->getMeshType();
      double zin = seg->getInletZ();
      double zout = seg->getOutletZ();

      cvOneDSubdomain* subdomain = new cvOneDSubdomain;
      assert(subdomain != 0);
      subdomain -> SetNumberOfNodes(nels+1);
      subdomain -> SetNumberOfElements(nels);
      subdomain -> SetMeshType(mType);
      subdomain -> Init(zin, zout);


      // Get the Initial Properties of the subdomain...
      double Qo = 0.0;
      double P0 = 0.0;
      double dQ0_dT = 0.0;
      if(i == 0) {
        Qo = seg->getInitialFlow();
        P0 =  seg->getInitialPressure();
        dQ0_dT = 0.0;
      }
      double So = seg->getInitInletS();
      double Sn = seg->getInitOutletS();
      BoundCondType boundT = seg -> getBoundCondition();
      double  boundV= seg -> getBoundValue();

      // Set these in the subdomain.
      subdomain->SetInitialFlow(Qo);
      subdomain->SetInitialdFlowdT(dQ0_dT);
      subdomain->SetInitInletS(So);
      subdomain->SetInitialPressure(P0);
      subdomain->SetInitOutletS(Sn);
      subdomain->SetGlobal1stNodeID(temp);
      subdomain->SetBoundCondition(boundT);
      if(!seg->IsOutlet){
        subdomain->SetBoundCondition(BoundCondTypeScope::NOBOUND);
      }
      subdomain->SetupMaterial(matID);

      // Set up Minor Loss
      subdomain->SetMinorLossType(seg->GetMinorLossType());
      if(seg->GetMinorLossType() != MinorLossScope::NONE){
        subdomain->SetBranchAngle(seg->GetBranchAngle());
        subdomain->SetUpstreamSeg(seg->GetUpstreamSeg());
        subdomain->SetBranchSeg(seg->GetBranchSeg());
      }

      // Set up boundary condition
      if (boundT == BoundCondTypeScope::RESISTANCE_TIME){
        double* time;
        double* resist;
        int num; // x-fer info from Segment to subdomain
        seg->getBoundResistanceValues(&resist,&time,&num);
        subdomain->SetBoundResistanceWave(time, resist, num);
      }else if(boundT == BoundCondTypeScope::RCR){
        double* rcr;
        int num;
        seg->getBoundRCRValues(&rcr,&num);
        subdomain -> SetBoundRCRValues(rcr,num);
        // cout<<"RCR boundary condition"<<endl;

      }else if(boundT == BoundCondTypeScope::RESISTANCE){ // modified resistance taking Pd wgyang
        double* resistance_pd;
        int num;
        seg->getBoundRCRValues(&resistance_pd,&num);
        subdomain -> SetBoundResistPdValues(resistance_pd,num);
        // cout<<"RESISTANCE boundary condition"<<endl;

      }else if (boundT == BoundCondTypeScope::PRESSURE_WAVE){
        double* time;
        double* pres;
        int num;
        seg->getBoundPressureValues(&pres,&time,&num);
        subdomain ->SetBoundPresWave(time, pres, num);

      }else if(boundT == BoundCondTypeScope::CORONARY){// Jongmin & Hyunjin
        double* time;
        double* p_lv;
        int num;
        seg->getBoundCoronaryValues(&p_lv, &time, &num);
        subdomain->SetBoundCoronaryValues(time, p_lv,num);
        // cout<<"CORONARY boundary condition"<<endl;

      }else{
        subdomain -> SetBoundValue(boundV);
      }


      temp += nels + 1;
      subdomainList.push_back(subdomain);
      if(seg->getIsOutletInfo()){
        outletList.push_back(i);
      }
      if(seg->getIsOutletInfo() == false){
        subdomain->SetBoundCondition(BoundCondTypeScope::NOBOUND);
      }
    }

    // For branch use.
    temp *= 2;
    for(i=0; i<ij; i++){
      jointList[i]->SetGlobal1stLagNodeID(temp);
      temp += jointList[i]->getNumberOfSegments();
    }
}

void cpl1DType::CreateGlobalArrays(void){
    assert( wasSet == false);
    long neq = mathModels[0]->GetTotalNumberOfEquations();

    long* maxa = new long[neq + 1];
    assert( maxa != 0);
    clear( neq + 1, maxa);

    int neeq = mathModels[0]->GetNumberOfElementEquations();
    // cout <<"Number of equations " << neq << endl;
    long* eqNumbers = new long[neeq];
    assert( eqNumbers != 0);
    long minEq, total;
    int i;

    // Element nodes
    total = 0;
    // std::cout << "subdomainList.size() = " << subdomainList.size() << std::endl;
    for(int i = 0; i < subdomainList.size(); i++){
      total += subdomainList[i]->GetNumberOfNodes();
      for( long el = 0; el < subdomainList[i]->GetNumberOfElements(); el++){
        mathModels[0]->GetEquationNumbers(el, eqNumbers, i);
        minEq = min(neeq, eqNumbers);
        for( int k = 0; k < neeq; k++){
          long currHeight = maxa[eqNumbers[k]];
          maxa[eqNumbers[k]] = max(currHeight, eqNumbers[k] - minEq);
        }
      }
    }

    // Joints
    total *= 2;
    for(i = 0; i < jointList.size(); i++){
      for(int j = 0; j < jointList[i]->getNumberOfSegments(); j++){
        minEq = mathModels[1]->GetUpmostEqnNumber(j, i);
        maxa[total] = total - minEq;
        total ++;
      }
    }

    // Now maxa contains the column heights
    // change it to hold the position
    // of the first element of the skyline at each column
    maxa[neq] = sum(neq, maxa);
    for( i = neq - 1; i >= 0; i--)
        maxa[i] = maxa[i+1] - maxa[i];

    // INITIALIZE MATRIX STORAGE SCHEME
    // AND ASSOCIATED SOLVER
# ifdef USE_SKYLINE
    lhs = new cvOneDSkylineMatrix(neq, maxa, "globalMatrix");
    // cvOneDGlobal::solver = new cvOneDSkylineLinearSolver();
    mathModels[0]->solver = new cvOneDSkylineLinearSolver();
# endif

# ifdef USE_SUPERLU
    lhs = new cvOneDSparseMatrix(neq, maxa, "globalMatrix");
    // cvOneDGlobal::solver = new cvOneDSparseLinearSolver();
    mathModels[0]->solver = new cvOneDSparseLinearSolver();
# endif

# ifdef USE_CSPARSE
    lhs = new cvOneDSparseMatrix(neq, maxa, "globalMatrix");
    // cvOneDGlobal::solver = new cvOneDSparseLinearSolver();
    mathModels[0]->solver = new cvOneDSparseLinearSolver();
# endif

    assert(lhs != 0);

    rhs = new cvOneDFEAVector( neq);
    assert(rhs != 0);

    // cvOneDGlobal::solver->SetLHS(lhs);
    // cvOneDGlobal::solver->SetRHS(rhs);
    // 只改了mathModels[0]，也就是cvOneDMthSegmentModel的
    mathModels[0]->solver->lhsMatrix = lhs;
    mathModels[0]->solver->rhsVector = rhs;


    previousSolution = new cvOneDFEAVector(neq, "previousSolution");
    assert(previousSolution != 0);
    previousSolution->Clear();
    currentSolution = new cvOneDFEAVector(neq, "currentSolution");
    assert(currentSolution != 0);
    currentSolution->Clear();
    increment = new cvOneDFEAVector(neq, "increment");
    assert(increment != 0);
    increment->Clear();
}

// Initialize the solution, that is, area as area input and flow rate as 0 except the inlet
void cpl1DType::CalcInitProps(long ID){
  double segLen = subdomainList[ID] -> GetLength();
  double Qo, dQ0dT;
  Qo = subdomainList[ID] -> GetInitialFlow();
  dQ0dT=0;

  double So = subdomainList[ID] -> GetInitInletS();
  double Sn = subdomainList[ID] -> GetInitOutletS();
  for( long node = 0; node < subdomainList[ID]->GetNumberOfNodes(); node++){
  double zn = subdomainList[ID]->GetNodalCoordinate( node);
  long eqNumbers[2];
  mathModels[0]->GetNodalEquationNumbers(node, eqNumbers, ID);

  // Linear Interpolation
  double zi = (zn - segLen)/(0.0-segLen);
  double Si = (zi*(So - Sn)) + Sn;
  (*previousSolution)[eqNumbers[0]] = Si;

    if(node == 0){
      (*previousSolution)[eqNumbers[1]] = Qo;
    }else{
      (*previousSolution)[eqNumbers[1]] = 0.0;
    }
  }
}

// =================
// GENERATE SOLUTION
// =================
void cpl1DType::GenerateSolution(){

  // // Print the formulation used
  // if(cpl1DType::useIV){
  //   cout << "Using Conservative Form ..." << endl;
  // }else{
  //   cout << "Using Advective Form ..." << endl;
  // }

  // Allocate the TotalSolution Array.
  TotalSolution.SetSize(Virtual_max_time_step +1, currentSolution -> GetDimension());
  // cout << "Total Solution is: " << Virtual_max_time_step  << " x ";
  // cout << currentSolution -> GetDimension() << endl;

  previousSolution->Rename( "step_0");
  *currentSolution = *previousSolution;

  double* tmp = previousSolution -> GetEntries();
  for (int j=0;j<previousSolution -> GetDimension(); j++){
    TotalSolution[0][j] = tmp[j];
  }

  // Initialize the Equations...
  for(int i = 0; i < mathModels.size(); i++){
    mathModels[i]->EquationInitialize(previousSolution, currentSolution);
  }

  //设置输出流及其精度
  std::ofstream outFile(cpl1DType::OutputFile, std::ios::app);
  std::cout << std::scientific << std::setprecision(3);
  outFile << std::scientific << std::setprecision(3);;

  cvOneDString String1( "step_");
  char String2[] = "99999";
  cvOneDString title;

  step = 1;
  int iter_total = 0;
  double toleration;
  clock_t tstart;
  clock_t tend;
  double time_consumed_total;
  
  tstart=clock();

  while(true){

    int numMath = mathModels.size();
    double currentTime = step * Virtual_time_step_size;
    increment->Clear();
    
    for(int i = 0; i < numMath; i++){ 
      mathModels[i]->TimeUpdate(currentTime, Virtual_time_step_size);
    }

    // Newton-Raphson Iterations...
    int iter = 0;
    double normf = 1.0;
    double norms = 1.0;
    clock_t tend_iter;
    clock_t tstart_iter;
    double timeConsumed;

    while(true){  //iter循环

      tstart_iter=clock();
        
      for(int i = 0; i < numMath; i++){
        mathModels[i]->FormNewton(lhs, rhs);
      }

      // PRINT RHS BEFORE BC APP
      if(cvOneDGlobal::debugMode){
        printf("(Debug) Printing LHS and RHS...\n");
        ofstream ofsRHS;
        ofstream ofsLHS;
        ofsRHS.open("rhs_1.txt");
        ofsLHS.open("lhs_1.txt");
        ofsRHS<<" --- RHS: Before ApplyBoundaryConditions" << endl;
        rhs->print(ofsRHS);
        ofsLHS<<" --- LHS: Before ApplyBoundaryConditions" << endl;
        lhs->print(ofsLHS);
        ofsRHS.close();
        ofsLHS.close();
        printf("ECCOLO\n");
        getchar();
      }

      mathModels[0]->ApplyBoundaryConditions();

      // PRINT RHS AFTER BC APP
      if(cvOneDGlobal::debugMode){
        cout<<" --- RHS: After Application of BC " << endl;
        rhs->print(cout);
      }

      // Do not evaluate residuals of lagrange eqns
      if(jointList.size() != 0){
        normf = rhs->Norm(L2_norm,1,2, jointList[0]->GetGlobal1stLagNodeID());
        norms = rhs->Norm(L2_norm,0,2, jointList[0]->GetGlobal1stLagNodeID());
      }else{
        normf = rhs->Norm(L2_norm,1,2);
        norms = rhs->Norm(L2_norm,0,2);
      }

      if (std::isnan(norms) || std::isnan(normf)) {
          throw cvException("Calculated a NaN for the residual.");
      }
      
      // Check Newton-Raphson Convergence
      if((currentTime != Virtual_time_step_size || (currentTime == Virtual_time_step_size && iter != 0)) && normf < cpl1DType::convergenceTolerance && norms < cpl1DType::convergenceTolerance){
        tend_iter=clock();
        timeConsumed = ((float)(tend_iter-tstart_iter))/CLOCKS_PER_SEC;
        outFile << "[" << outletName << "]: " << step << "-" << iter << "  " << normf << "  " << norms << "  " << timeConsumed << endl;
        break;
      }

      // Add increment
      increment->Clear();

      // cvOneDGlobal::solver->Solve(*increment);
      mathModels[0]->solver->Solve(*increment);

      currentSolution->Add(*increment);

      // If the area goes less than zero, it tells in which segment the error occurs.
      // Assumes that all the lagrange multipliers are at the end of the vector.
      int negArea=0;
      if(jointList.size() != 0){
        for (long i= 0; i< jointList[0]->GetGlobal1stLagNodeID();i+=2){
          long elCount = 0;
          int fileIter = 0;
          //check if area <0 or =nan
          if (currentSolution->Get(i) < 0.0 || (currentSolution->Get(i) != currentSolution->Get(i))){
            negArea=1;
            while (fileIter < model -> getNumberOfSegments()){
              cvOneDSegment *curSeg = model -> getSegment(fileIter);
              long numEls = curSeg -> getNumElements();
              long startOut = elCount;
              long finishOut = elCount + ((numEls+1)*2);
              char *modelname;
              char *segname;
              if (startOut <= i && i <= finishOut) {
                  modelname = model-> getModelName();
                  segname = curSeg -> getSegmentName();
                  std::string msg = "ERROR: The area of segment '" + std::string(segname) + "' is negative.";
                  throw cvException(msg.c_str());
              }
              elCount += 2*(numEls+1);
              fileIter++;
              }
            }
          }
        }

      if(negArea==1) {
      std::string emptyString;
      postprocess_Text(emptyString);
      assert(0);
      }

      if(cvOneDGlobal::debugMode){
        printf("(Debug) Printing Solution...\n");
        ofstream ofs("solution.txt");
        for(int loopA=0;loopA<currentSolution->GetDimension();loopA++){
          ofs << to_string(loopA) << " " << currentSolution->Get(numMath) << endl;
        }
        ofs.close();
        getchar();
      }


      // A flag in case the cross sectional area is negative, but don't want to include the lagrange multipliers
      // Assumes that all the lagrange multipliers are at the end of the vector
      if(jointList.size() != 0){
        currentSolution->CheckPositive(0,2,jointList[0]->GetGlobal1stLagNodeID());
      }else{
        currentSolution->CheckPositive(0,2,currentSolution->GetDimension());
      }

      // Set Boundary Conditions
      // cvOneDMthModelBase::CurrentInletFlow = flowEachTime;
      mathModels[0]->CurrentInletFlow = flowEachTime;
      mathModels[0]->SetBoundaryConditions();

      if(iter > MAX_NONLINEAR_ITERATIONS){
        outFile << "Error: Newton not converged, exceed max iterations" << endl;
        outFile << "norm of Flow rate:" << normf << ", norm of Area:" << norms << endl;
        outFile << "----------------------------------------------------" << endl;
        break;
      }

      iter++;
    }// End linear while

    sprintf( String2, "%ld", (unsigned long)step);
    title = String1 + String2;
    currentSolution->Rename(title.data());

    double * tmp = currentSolution -> GetEntries();
    for(int j=0;j<currentSolution -> GetDimension(); j++){
      TotalSolution[step][j] = tmp[j];
    }

    *previousSolution = *currentSolution;
    iter_total += iter;

    if(step > 1){
      long dimension = currentSolution -> GetDimension();
      cvOneDMaterial* curMat = subdomainList[0]->GetMaterial();
      double preCurrent = (double) curMat->GetPressure((double)TotalSolution[step][0],0); //第0个点最后一个时间步的面积
      double preBefore = (double) curMat->GetPressure((double)TotalSolution[step-1][0],0);
      // cout << preCurrent << endl;
      // cout << preBefore << endl;
      toleration = std::abs((preCurrent - preBefore)/preBefore);
      
      if( toleration < Time_Pressure_Relative_Tolerations || step == Virtual_max_time_step){
        tend=clock();
        time_consumed_total = ((float)(tend-tstart))/CLOCKS_PER_SEC;
        cout << "[" << outletName << "]: " << step << "-" << iter_total << "  " << toleration << "  " << time_consumed_total << endl;
        outFile << "--------------------------------------------------------------" << endl;
        preFrom1DEachTime = preCurrent;  //输出第一个点的压强给三维
        break;
      }
    }
    step ++;
  
  } // End nonlinear loop

  outFile.close();
}


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

#include "cpl1DType.h"

using namespace std;

// ====================================
// GET DATA TABLE ENTRY FROM STRING KEY
// ====================================
int cpl1DType::getDataTableIDFromStringKey(string key){
  bool found = false;
  int count = 0;
  while((!found)&&(count<cvOneDGlobal::gDataTables.size())){
    found = upper_string(key) == upper_string(cvOneDGlobal::gDataTables[count]->getName());
    // Update Counter
    if(!found){
      count++;
    }
  }
  if(!found){
    throw cvException(string("ERROR: Cannot find data table entry from string key: " + key + ".\n").c_str());
    return -1;
  }else{
    return count;
  }
}

// ===============================
// CREATE MODEL AND RUN SIMULATION
// ===============================
void cpl1DType::createModel(){

  // CREATE MODEL MANAGER
  cvOneDModelManager* oned = new cvOneDModelManager((char*)opts.modelName.c_str());

  // CREATE NODES
  // cout << "Creating Nodes ... " << endl;
  int totNodes = opts.nodeName.size();
  int nodeError = CV_OK;
  for(int loopA=0;loopA<totNodes;loopA++){
    // Finally Create Joint
    nodeError = oned->CreateNode((char*)opts.nodeName[loopA].c_str(),
                                   opts.nodeXcoord[loopA], opts.nodeYcoord[loopA], opts.nodeZcoord[loopA]);
    if(nodeError == CV_ERROR){
      throw cvException(string("ERROR: Error Creating NODE " + to_string(loopA) + "\n").c_str());
    }
  }

  // CREATE JOINTS
  // cout << "Creating Joints ... " << endl;
  int totJoints = opts.jointName.size();
  int jointError = CV_OK;
  int* asInlets = NULL;
  int* asOutlets = NULL;
  string currInletName;
  string currOutletName;
  int jointInletID = 0;
  int jointOutletID = 0;
  int totJointInlets = 0;
  int totJointOutlets = 0;
  for(int loopA=0;loopA<totJoints;loopA++){
    // GET NAMES FOR INLET AND OUTLET
    currInletName = opts.jointInletName[loopA];
    currOutletName = opts.jointOutletName[loopA];
    // FIND JOINTINLET INDEX
    jointInletID = getListIDWithStringKey(currInletName,opts.jointInletListNames);
    if(jointInletID < 0){
      throw cvException(string("ERROR: Cannot Find JOINTINLET for key " + currInletName).c_str());
    }
    totJointInlets = opts.jointInletListNumber[jointInletID];
    // FIND JOINTOUTLET INDEX
    jointOutletID = getListIDWithStringKey(currOutletName,opts.jointOutletListNames);
    if(jointInletID < 0){
      throw cvException(string("ERROR: Cannot Find JOINTOUTLET for key " + currOutletName).c_str());
    }
    // GET TOTALS
    totJointInlets = opts.jointInletListNumber[jointInletID];
    totJointOutlets = opts.jointOutletListNumber[jointOutletID];
    // ALLOCATE INLETS AND OUTLET LIST
    asInlets = NULL;
    asOutlets = NULL;
    if(totJointInlets > 0){
      asInlets = new int[totJointInlets];
      for(int loopB=0;loopB<totJointInlets;loopB++){
        asInlets[loopB] = opts.jointInletList[jointInletID][loopB];
      }
    }
    if(totJointOutlets > 0){
      asOutlets = new int[totJointOutlets];
      for(int loopB=0;loopB<totJointOutlets;loopB++){
        asOutlets[loopB] = opts.jointOutletList[jointOutletID][loopB];
      }
    }
    // Finally Create Joint
    jointError = oned->CreateJoint((char*)opts.jointName[loopA].c_str(),
                                   opts.nodeXcoord[loopA], opts.nodeYcoord[loopA], opts.nodeZcoord[loopA],
                                   totJointInlets, totJointOutlets,asInlets,asOutlets);
    if(jointError == CV_ERROR){
      throw cvException(string("ERROR: Error Creating JOINT " + to_string(loopA) + "\n").c_str());
    }
    // Deallocate
    delete [] asInlets;
    delete [] asOutlets;
    asInlets = NULL;
    asOutlets = NULL;
  }

  // CREATE MATERIAL
  // cout << "Creating Materials ... " << endl;
  int totMaterials = cpl1DType::materialName.size();
  // cout << "totMaterials = " << totMaterials << endl;
  //在每个处理cpl1D的进程上都创建所有的MATERIAL以便使用
  int matError = CV_OK;
  double doubleParams[3];
  int matID = 0;
  string currMatType = "MATERIAL";
  int numParams = 0;
  for(int loopA = 0;loopA<totMaterials;loopA++){
    if(upper_string(cpl1DType::materialType[loopA]) == "OLUFSEN"){
      currMatType = "MATERIAL_OLUFSEN";
      numParams = 3;
    }else{
      currMatType = "MATERIAL_LINEAR";
      numParams = 1;
    }
    doubleParams[0] = cpl1DType::materialParam1[loopA];
    doubleParams[1] = cpl1DType::materialParam2[loopA];
    doubleParams[2] = cpl1DType::materialParam3[loopA];
    // CREATE MATERIAL
    matError = oned->CreateMaterial((char*)cpl1DType::materialName[loopA].c_str(),
                                    (char*)currMatType.c_str(),
                                    cpl1DType::materialDensity[loopA],
                                    cpl1DType::materialViscosity[loopA],
                                    cpl1DType::materialExponent[loopA],
                                    cpl1DType::materialPRef[loopA],
                                    numParams, doubleParams,
                                    &matID);
    if(matError == CV_ERROR){
      throw cvException(string("ERROR: Error Creating MATERIAL " + to_string(loopA) + "\n").c_str());
    }

  }

  // CREATE DATATABLES
  // cout << "Creating Data Tables ... " << endl;
  int totCurves = opts.dataTableName.size();
  int curveError = CV_OK;
  for(int loopA=0;loopA<totCurves;loopA++){
    curveError = oned->CreateDataTable((char*)opts.dataTableName[loopA].c_str(),(char*)opts.dataTableType[loopA].c_str(),opts.dataTableVals[loopA]);
    if(curveError == CV_ERROR){
      throw cvException(string("ERROR: Error Creating DATATABLE " + to_string(loopA) + "\n").c_str());
    }
  }

  // SEGMENT DATA
  // cout << "Creating Segments ... " << endl;
  int segmentError = CV_OK;
  int totalSegments = opts.segmentName.size();
  int curveTotals = 0;
  double* curveTime = NULL;
  double* curveValue = NULL;
  string matName;
  string curveName;
  int currMatID = 0;
  int dtIDX = 0;
  for(int loopA=0;loopA<totalSegments;loopA++){

    // GET MATERIAL
    matName = opts.segmentMatName[loopA];
    currMatID = getListIDWithStringKey(matName,cpl1DType::materialName);
    // cout << "matName = " << matName << endl;
    if(currMatID < 0){
      throw cvException(string("ERROR: Cannot Find Material for key " + matName).c_str());
    }

    // GET CURVE DATA
    curveName = opts.segmentDataTableName[loopA];

    if(upper_string(curveName) != "NONE") {
      dtIDX = getDataTableIDFromStringKey(curveName);
      curveTotals = cvOneDGlobal::gDataTables[dtIDX]->getSize();
      curveTime = new double[curveTotals];
      curveValue = new double[curveTotals];
      for(int loopA=0;loopA<curveTotals;loopA++){
        curveTime[loopA] = cvOneDGlobal::gDataTables[dtIDX]->getTime(loopA);
        curveValue[loopA] = cvOneDGlobal::gDataTables[dtIDX]->getValues(loopA);
      }
    }else{
      curveTotals = 1;
      curveTime = new double[curveTotals];
      curveValue = new double[curveTotals];
      curveTime[0] = 0.0;
      curveValue[0] = 0.0;
    }
    segmentError = oned->CreateSegment((char*)opts.segmentName[loopA].c_str(),
                                       (long)opts.segmentID[loopA],
                                       opts.segmentLength[loopA],
                                       (long)opts.segmentTotEls[loopA],
                                       (long)opts.segmentInNode[loopA],
                                       (long)opts.segmentOutNode[loopA],
                                       opts.segmentInInletArea[loopA],
                                       opts.segmentInOutletArea[loopA],
                                       opts.segmentInFlow[loopA],
                                       currMatID,
                                       (char*)opts.segmentLossType[loopA].c_str(),
                                       opts.segmentBranchAngle[loopA],
                                       opts.segmentUpstreamSegment[loopA],
                                       opts.segmentBranchSegment[loopA],
                                       (char*)opts.segmentBoundType[loopA].c_str(),
                                       curveValue,
                                       curveTime,
                                       curveTotals);
    if(segmentError == CV_ERROR){
      throw cvException(string("ERROR: Error Creating SEGMENT " + to_string(loopA) + "\n").c_str());
    }
    // Deallocate
    delete [] curveTime;
    curveTime = NULL;
    delete [] curveValue;
    curveValue = NULL;
  }

  model = cvOneDGlobal::gModelList[cvOneDGlobal::currentModel];
  prepro();

}

// ======================
// READ SINGLE MODEL FILE
// ======================
void cpl1DType::readModelFile(cvStringVec includedFiles){

  // Message
  // cout << endl;
  // cout << "Reading file: " << inputFile.c_str() << "..." << endl;

  // Declare
  cvStringVec tokenizedString;
  cvLongVec   tempIntVec;
  string      matType;
  cvDoubleVec temp;
  bool doInclude = false;

  // Declare input File
  ifstream infile;
  infile.open(inputFile);
  if(infile.fail()){
    throw cvException("ERROR: Input file does not exist.\n");
  }
  int lineCount = 1;
  int reminder = 0;
  int totSegments = 0;

  // Read Data From File
  std::string buffer;
  while (std::getline(infile,buffer)){

    // Trim String
    buffer = trim_string(buffer);

    // Tokenize String
    tokenizedString = split_string(buffer, " ,\t");
    if (tokenizedString.size() == 0) { 
      continue;
    }

    // Check for Empty buffer
    if(!buffer.empty()){
      // CHECK THE ELEMENT TYPE
      if(upper_string(tokenizedString[0]) == "MODEL"){
        //cout << "Found Model.\n");
        if(opts.modelNameDefined){
          throw cvException("ERROR: Model name already defined\n");
        }
        if(tokenizedString.size() > 2){
          throw cvException(string("ERROR: Too many parameters for MODEL token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 2){
          throw cvException(string("ERROR: Not enough parameters for MODEL token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          opts.modelName = tokenizedString[1];
        }catch(...){
          throw cvException(string("ERROR: Invalid Model Name. Line " + to_string(lineCount) + "\n").c_str());
        }
        opts.modelNameDefined = true;
      }else if(upper_string(tokenizedString[0]) == "NODE"){
        // cout << "Found Joint.\n");
        if(tokenizedString.size() > 5){
          throw cvException(string("ERROR: Too many parameters for NODE token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 5){
          throw cvException(string("ERROR: Not enough parameters for NODE token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // Get Node Name
          opts.nodeName.push_back(tokenizedString[1]);
          opts.nodeXcoord.push_back(atof(tokenizedString[2].c_str()));
          opts.nodeYcoord.push_back(atof(tokenizedString[3].c_str()));
          opts.nodeZcoord.push_back(atof(tokenizedString[4].c_str()));
        }catch(...){
          throw cvException(string("ERROR: Invalid NODE Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(upper_string(tokenizedString[0]) == "JOINT"){
        // cout << "Found Joint.\n");
        if(tokenizedString.size() > 5){
          throw cvException(string("ERROR: Too many parameters for JOINT token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 5){
          throw cvException(string("ERROR: Not enough parameters for JOINT token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // Get Joint Name
          opts.jointName.push_back(tokenizedString[1]);
          //opts.jointNode.push_back(atof(tokenizedString[2].c_str()));
          opts.jointNode.push_back(tokenizedString[2]);
          opts.jointInletName.push_back(tokenizedString[3]);
          opts.jointOutletName.push_back(tokenizedString[4]);
        }catch(...){
          throw cvException(string("ERROR: Invalid JOINT Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(upper_string(tokenizedString[0]) == "JOINTINLET"){
        // cout << "Found JointInlet.\n");
        if(tokenizedString.size() < 3){
          throw cvException(string("ERROR: Not enough parameters for JOINTINLET token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // Get Joint Name
          opts.jointInletListNames.push_back(tokenizedString[1]);
          totSegments = atoi(tokenizedString[2].c_str());
          opts.jointInletListNumber.push_back(totSegments);
          tempIntVec.clear();
          if(totSegments > 0){
            for(size_t loopA=3;loopA<tokenizedString.size();loopA++){
              tempIntVec.push_back(atoi(tokenizedString[loopA].c_str()));
            }
          }
          opts.jointInletList.push_back(tempIntVec);
        }catch(...){
          throw cvException(string("ERROR: Invalid JOINTINLET Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(upper_string(tokenizedString[0]) == "JOINTOUTLET"){
        // cout << "Found JointOutlet.\n");
        if(tokenizedString.size() < 3){
          throw cvException(string("ERROR: Not enough parameters for JOINTOUTLET token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // Get Joint Name
          opts.jointOutletListNames.push_back(tokenizedString[1]);
          totSegments = atoi(tokenizedString[2].c_str());
          opts.jointOutletListNumber.push_back(totSegments);
          tempIntVec.clear();
          if(totSegments > 0){
            for(size_t loopA=3;loopA<tokenizedString.size();loopA++){
              tempIntVec.push_back(atoi(tokenizedString[loopA].c_str()));
            }
          }
          opts.jointOutletList.push_back(tempIntVec);
        }catch(...){
          throw cvException(string("ERROR: Invalid JOINTOUTLET Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(upper_string(tokenizedString[0]) == "SEGMENT"){
        // cout << "Found Segment.\n");
        if(tokenizedString.size() > 17){
          throw cvException(string("ERROR: Too many parameters for SEGMENT token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 17){
          throw cvException(string("ERROR: Not enough parameters for SEGMENT token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // char* segName,
          opts.segmentName.push_back(tokenizedString[1]);
          // long segID,
          opts.segmentID.push_back(atoi(tokenizedString[2].c_str()));
          // double  segLen,
          opts.segmentLength.push_back(atof(tokenizedString[3].c_str()));
          // long    numEls,
          opts.segmentTotEls.push_back(atoi(tokenizedString[4].c_str()));
          // long    inNode,
          opts.segmentInNode.push_back(atoi(tokenizedString[5].c_str()));
          // long    outNode,
          opts.segmentOutNode.push_back(atoi(tokenizedString[6].c_str()));
          // double  InitialInletArea,
          opts.segmentInInletArea.push_back(atof(tokenizedString[7].c_str()));
          // double  InitialOutletArea,
          opts.segmentInOutletArea.push_back(atof(tokenizedString[8].c_str()));
          // double  InitialFlow,
          opts.segmentInFlow.push_back(atof(tokenizedString[9].c_str()));
          // int matID,
          opts.segmentMatName.push_back(tokenizedString[10].c_str());
          // char* lossType,
          opts.segmentLossType.push_back(tokenizedString[11]);
          // double branchAngle,
          opts.segmentBranchAngle.push_back(atof(tokenizedString[12].c_str()));
          // int upstreamSegment,
          opts.segmentUpstreamSegment.push_back(atoi(tokenizedString[13].c_str()));
          // int branchSegment,
          opts.segmentBranchSegment.push_back(atoi(tokenizedString[14].c_str()));
          // char* boundType,
          opts.segmentBoundType.push_back(tokenizedString[15]);
          // Curve ID Instead of num,value,time
          // double* value,
          // double* time,
          // int num
          opts.segmentDataTableName.push_back(tokenizedString[16]);
        }catch(...){
          throw cvException(string("ERROR: Invalid SEGMENT Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(upper_string(tokenizedString[0]) == "DATATABLE"){
        // cout << "Found Data Table.\n");
        try{

          // Get Datatable Name
          opts.dataTableName.push_back(tokenizedString[1]);
          // Add the type of the datatable
          opts.dataTableType.push_back(tokenizedString[2]);

          bool foundEnd = false;
          temp.clear();
          while(!foundEnd){
            std::getline(infile,buffer);
            lineCount++;
            // Trim String
            buffer = trim_string(buffer);
            // Tokenize String
            tokenizedString = split_string(buffer, " ,\t");
            if (tokenizedString.size() == 0) { 
              break;
            }
            // Check for Empty buffer
            if(!buffer.empty()){
              if(upper_string(tokenizedString[0]) == std::string("ENDDATATABLE")){
                foundEnd = true;
              }else{
                for(int loopA=0;loopA<tokenizedString.size();loopA++){
                  temp.push_back(atof(tokenizedString[loopA].c_str()));
                }
              }
            }
          }
          // Add all the values to the option array
          opts.dataTableVals.push_back(temp);
        }catch(...){
          throw cvException(string("ERROR: Invalid DATATABLE Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(upper_string(tokenizedString[0]) == "INCLUDE"){
        // Check if the file is active
        if(upper_string(tokenizedString[2]) == "TRUE"){
          doInclude = true;
        }else if(upper_string(tokenizedString[2]) == "FALSE"){
          doInclude = false;
        }else{
          throw cvException(string("ERROR: Invalid INCLUDE switch format. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // If active include file in list
          if(doInclude){
            includedFiles.push_back(tokenizedString[1]);
          }
        }catch(...){
          throw cvException(string("ERROR: Invalid INCLUDE Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if((tokenizedString.size() == 0)||(tokenizedString[0].at(0) == '#')||(tokenizedString[0].find_first_not_of(' ') == std::string::npos)){
        // cout << "Found Blank.\n");
        // COMMENT OR BLANK LINE: DO NOTHING
      }else if(upper_string(tokenizedString[0]) == "SOLVEROPTIONS"){
      }else if(upper_string(tokenizedString[0]) == "MATERIAL"){
      }else if(upper_string(tokenizedString[0]) == "OUTPUT"){
      }else{
        // TOKEN NOT RECOGNIZED
        throw cvException(string("ERROR: Invalid Token in input file, line: "  + to_string(lineCount) + "\n").c_str());
      }
    }
    // cout << "Line: %d, Buffer: %s\n",lineCount,buffer.c_str());
    // getchar();

    // Increment Line Count
    lineCount++;
  }
  // Close File
  infile.close();
}

// ====================
// READ MODEL FROM FILE
// ====================
void cpl1DType::readModel(){

  // List of included Files
  cvStringVec includedFiles;
  string currentFile;

  // Read First File
  readModelFile(includedFiles);

  //Read Nested Files
  while(includedFiles.size() > 0){

    // Get the first file Name
    currentFile = includedFiles[0];
    // Delete the First element
    includedFiles.erase(includedFiles.begin());
    // Read the file and store new included files
    readModelFile(includedFiles);
  }
}

void cpl1DType::readSharedVar(){
  // Declare input File
  ifstream infile;
  infile.open(inputFile);
  if(infile.fail()){
    throw cvException("ERROR: Input file does not exist.\n");
  }

  // Read Data From File
  std::string buffer;
  cvStringVec tokenizedString;
  int lineCount = 1;
  string  matType;
  while (std::getline(infile,buffer)){

    // Trim String
    buffer = trim_string(buffer);

    // Tokenize String
    tokenizedString = split_string(buffer, " ,\t");
    if (tokenizedString.size() == 0) { 
      continue;
    }
    // Check for Empty buffer
    if(!buffer.empty()){
      // CHECK THE ELEMENT TYPE
      if(upper_string(tokenizedString[0]) == "SOLVEROPTIONS"){
        // cout << "Found Solver Options." << endl;
        if(solverOptionDefined){
          throw cvException("ERROR: SOLVEROPTIONS already defined\n");
        }
        if(tokenizedString.size() > 5){
          throw cvException(string("ERROR: Too many parameters for SOLVEROPTIONS token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 5){
          throw cvException(string("ERROR: Not enough parameters for SOLVEROPTIONS token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // long quadPoints,
          quadPoints = atoi(tokenizedString[1].c_str());
          // double conv,
          convergenceTolerance = atof(tokenizedString[2].c_str());
          // int useIV,
          useIV = atoi(tokenizedString[3].c_str());
          // int usestab
          useStab = atoi(tokenizedString[4].c_str());
        }catch(...){
          throw cvException(string("ERROR: Invalid SOLVEROPTIONS Format. Line " + to_string(lineCount) + "\n").c_str());
        }
        solverOptionDefined = true;
      }else if(upper_string(tokenizedString[0]) == "OUTPUT"){
        if(tokenizedString.size() > 2){
          throw cvException(string("ERROR: Too many parameters for OUTPUT token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 2){
          throw cvException(string("ERROR: Not enough parameters for OUTPUT token. Line " + to_string(lineCount) + "\n").c_str());
        }
        // Output Type
        if(upper_string(tokenizedString[1]) == "TEXT"){
          outputType = OutputTypeScope::OUTPUT_TEXT;
        }else if(upper_string(tokenizedString[1]) == "VTK"){
          outputType = OutputTypeScope::OUTPUT_VTK;
        }else if(upper_string(tokenizedString[1]) == "BOTH"){
          outputType = OutputTypeScope::OUTPUT_BOTH;
        }else{
          throw cvException("ERROR: Invalid OUTPUT Type.\n");
        }
      }else if(upper_string(tokenizedString[0]) == "MATERIAL"){
        // cout << "Found Material." << endl;
        opts.MATnum += 1;
        if(tokenizedString.size() > 10){
          throw cvException(string("ERROR: Too many parameters for MATERIAL token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 8){
          throw cvException(string("ERROR: Not enough parameters for MATERIAL token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // Material Name
          materialName.push_back(tokenizedString[1]);
          // Material Type
          matType = tokenizedString[2];
          materialType.push_back(matType);
          // Density
          materialDensity.push_back(atof(tokenizedString[3].c_str()));
          // Dynamic Viscosity
          materialViscosity.push_back(atof(tokenizedString[4].c_str()));
          // Reference Pressure
          materialPRef.push_back(atof(tokenizedString[5].c_str()));
          // Material Exponent
          materialExponent.push_back(atof(tokenizedString[6].c_str()));
          // Extra Material Parameters
          if(upper_string(matType) == "OLUFSEN"){
            materialParam1.push_back(atof(tokenizedString[7].c_str()));
            materialParam2.push_back(atof(tokenizedString[8].c_str()));
            materialParam3.push_back(atof(tokenizedString[9].c_str()));
          }else if(upper_string(matType) == "LINEAR"){
            materialParam1.push_back(atof(tokenizedString[7].c_str()));
            materialParam2.push_back(0.0);
            materialParam3.push_back(0.0);
          }else{
            throw cvException(string("ERROR: Invalid MATERIAL Type. Line " + to_string(lineCount) + "\n").c_str());
          }
        }catch(...){
          throw cvException("ERROR: Invalid MATERIAL Format.\n");
        }
      }else if(upper_string(tokenizedString[0]) == "MODEL"){
      }else if(upper_string(tokenizedString[0]) == "NODE"){
      }else if(upper_string(tokenizedString[0]) == "JOINT"){
      }else if(upper_string(tokenizedString[0]) == "JOINTINLET"){
      }else if(upper_string(tokenizedString[0]) == "JOINTOUTLET"){
      }else if(upper_string(tokenizedString[0]) == "SEGMENT"){
      }else if(upper_string(tokenizedString[0]) == "DATATABLE"){
        bool foundEnd = false;
          while(!foundEnd){
            std::getline(infile,buffer);
            lineCount++;
            // Trim String
            buffer = trim_string(buffer);
            // Tokenize String
            tokenizedString = split_string(buffer, " ,\t");
            if (tokenizedString.size() == 0) { 
              break;
            }
            // Check for Empty buffer
            if(!buffer.empty()){
              if(upper_string(tokenizedString[0]) == std::string("ENDDATATABLE")){
                foundEnd = true;
              }
            }
          }
      }else if(upper_string(tokenizedString[0]) == "ENDDATATABLE"){
      }else if(upper_string(tokenizedString[0]) == "INCLUDE"){
      }else if((tokenizedString.size() == 0)||(tokenizedString[0].at(0) == '#')||(tokenizedString[0].find_first_not_of(' ') == std::string::npos)){
      }else{
        // COMMENT OR BLANK LINE: DO NOTHING
        throw cvException(string("ERROR: Invalid Token in input file, line: "  + to_string(lineCount) + "\n").c_str());
      }
    }
    lineCount++;
  }
  // Close File
  infile.close();
  // Check for double material name
  int dblIDX = -1;
  // PERFORM DATA CHECKING
  for(int loopA=0;loopA < materialName.size();loopA++){
    for(int loopB=loopA+1;loopB < materialName.size();loopB++){
      if(materialName[loopA] == materialName[loopB]){
        dblIDX = loopA;
      }
    }
  }
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double Material Name: " + materialName[dblIDX] + "\n").c_str());
  }
}
