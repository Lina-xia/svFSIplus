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

#include "cvOneDOptions.h"

double cvOneDOptions::dt = 0;
int   cvOneDOptions::saveIncr = 0;  //incre
int   cvOneDOptions::maxStep = 0;
int   cvOneDOptions::quadPoints = 2;
double cvOneDOptions::convergenceTolerance = 1.0e-8;
int    cvOneDOptions::useIV = 1;
int    cvOneDOptions::useStab = 1;
int    cvOneDOptions::outputType = 0; // Default Text Output
int    cvOneDOptions::vtkOutputType = 0; // Default Multiple Files

// CONSTRUCTOR
cvOneDOptions::cvOneDOptions(){
  modelNameDefined = false;
  solverOptionDefined = false;
}

// DESTRUCTOR
cvOneDOptions::~cvOneDOptions(){
}

void cvOneDOptions::checkSegmentLengthConsistency(){
  bool inconsistencyFound = false;
  // Get total number of segments
  int totSegs = segmentLength.size();
  double segLength = 0.0;
  int inNode = 0;
  int outNode = 0;
  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
  double nodeDist = 0.0;
  for(int loopA=0;loopA<totSegs;loopA++){
    // Get Current Segment Length
    segLength = segmentLength[loopA];
    // Get end nodes
    inNode = segmentInNode[loopA];
    outNode = segmentOutNode[loopA];
    // Get node spatial distance
    dx = nodeXcoord[outNode] - nodeXcoord[inNode];
    dy = nodeYcoord[outNode] - nodeYcoord[inNode];
    dz = nodeZcoord[outNode] - nodeZcoord[inNode];

    nodeDist = sqrt(dx*dx + dy*dy + dz*dz);
  }
  if(inconsistencyFound){
    printf("WARNING: Inconsistency detected between segment length and end node distance.\n");
    printf("Changing the segment lengths.\n");
  }
}


// PERFORM DATA CHECKING
template <class T>
int checkForDoubleEntry(vector<T> vec){
  for(int loopA=0;loopA<vec.size();loopA++){
    for(int loopB=loopA+1;loopB<vec.size();loopB++){
      if(vec[loopA] == vec[loopB]){
        return loopA;
      }
    }
  }
  return -1;
}

// PERFORM DATA CHECKING
template <class T>
bool checkContains(T value, vector<T> vec){
  for(int loopA=0;loopA<vec.size();loopA++){
    if(vec[loopA] == value){
      return true;
    }
  }
  return false;
}

// CHECK FOR POSITIVE VALUES
int checkForPositiveVal(cvDoubleVec vec){
  for(int loopA=0;loopA<vec.size();loopA++){
    if(vec[loopA] < 0.0){
      return loopA;
    }
  }
  return -1;
}

// =====================
// PERFORM DATA CHECKING
// =====================
void cvOneDOptions::check(){

  int dblIDX = 0;
  int chkIDX = 0;

  // Check for double node name
  dblIDX = checkForDoubleEntry(nodeName);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double Node Name: " + nodeName[dblIDX] + "\n").c_str());
  }
  // Check for double joint name
  dblIDX = checkForDoubleEntry(jointName);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double Joint Name: " + jointName[dblIDX] + "\n").c_str());
  }
  // Check for double jointInlet name
  dblIDX = checkForDoubleEntry(jointInletListNames);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double JointInlet Name: " + jointInletListNames[dblIDX] + "\n").c_str());
  }
   // Check for double jointOutlet name
  dblIDX = checkForDoubleEntry(jointOutletListNames);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double JointOutlet Name: " + jointOutletListNames[dblIDX] + "\n").c_str());
  }
  // Check for double material name
  dblIDX = checkForDoubleEntry(materialName);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double Material Name: " + materialName[dblIDX] + "\n").c_str());
  }
  // Check for double data table name
  dblIDX = checkForDoubleEntry(dataTableName);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double Data Table Name: " + dataTableName[dblIDX] + "\n").c_str());
  }
  // Check for double segment Name
  dblIDX = checkForDoubleEntry(segmentName);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double Segment Name: " + segmentName[dblIDX] + "\n").c_str());
  }
  // Check for double segment ID
  dblIDX = checkForDoubleEntry(segmentID);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double Segment ID: " + to_string(segmentID[dblIDX]) + "\n").c_str());
  }
  // Check for negative area in input
  chkIDX = checkForPositiveVal(segmentInInletArea);
  if(chkIDX >= 0){
    throw cvException(string("ERROR: Negative Inlet area in segment : " + segmentName[chkIDX] + "\n").c_str());
  }
  chkIDX = checkForPositiveVal(segmentInOutletArea);
  if(chkIDX >= 0){
    throw cvException(string("ERROR: Negative Outlet area in segment : " + segmentName[chkIDX] + "\n").c_str());
  }

  // Check the consistency between node coords and segment lengths
  checkSegmentLengthConsistency();

  // Check if the segments refer to a node that is not there
  checkSegmentHasNodes();

  // Check if the joints refer to a node that is not there
  checkJointHasNodes();

}

void cvOneDOptions::checkSegmentHasNodes(){
  int inNode = 0;
  int outNode = 0;
  for(int loopA=0;loopA<segmentName.size();loopA++){
    // Get end nodes
    inNode = segmentInNode[loopA];
    outNode = segmentOutNode[loopA];
    // Check
    if((inNode < 0)||(inNode >= nodeName.size())){
      throw cvException(string("ERROR: Missing Node in Segment: " + segmentName[loopA] + "\n").c_str());
    }
    if((outNode < 0)||(outNode >= nodeName.size())){
      throw cvException(string("ERROR: Missing Node in Segment: " + segmentName[loopA] + "\n").c_str());
    }
  }
}

void cvOneDOptions::checkJointHasNodes(){
  string currNodeName;
  for(int loopA=0;loopA<jointName.size();loopA++){
    // Get end nodes
    currNodeName = jointNode[loopA];
    // Check If
    if(!checkContains(currNodeName,nodeName)){
      throw cvException(string("ERROR: Missing Node " + currNodeName + " in Joint: " + jointName[loopA] + "\n").c_str());
    }
  }
}
