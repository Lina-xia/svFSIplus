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

# include "cvOneDEnums.h"
# include "cvOneDGlobal.h"
# include "cvOneDMaterialManager.h"
# include "cvOneDMaterialOlufsen.h"
# include "cvOneDMaterialLinear.h"


// Constructor
cvOneDMaterialManager::cvOneDMaterialManager(){
  numMaterials = 0;
  for (int i=0;i<MAX_NUM_OF_MATERIAL_INSTANCES;i++){
    materials[i] = NULL;
  }
}

// Destructor
cvOneDMaterialManager::~cvOneDMaterialManager(){
}

int cvOneDMaterialManager::AddNewMaterial(MaterialType type, cvOneDMaterial* mat){
  int rtnID = numMaterials;
  types[rtnID] = type;
  materials[rtnID] = mat;
  numMaterials++;
  return rtnID;
}

int cvOneDMaterialManager::AddNewMaterialOlufsen(double density, double dynamicViscosity,
                                                 double profile_exponent, double pRef,
                                                 double *params){
  cvOneDMaterialOlufsen* olfmat = new cvOneDMaterialOlufsen();
  olfmat->SetDensity(density);
  olfmat->SetDynamicViscosity(dynamicViscosity);
  olfmat->SetProfileExponent(profile_exponent);
  olfmat->SetReferencePressure(pRef);
  olfmat->SetMaterialType(params,pRef);
  // cout << "new cvOneMaterialOlufsen called check pRef " << olfmat->GetReferencePressure() << endl;
  return cvOneDGlobal::gMaterialManager->AddNewMaterial(MaterialType_MATERIAL_OLUFSEN,(cvOneDMaterial*)olfmat);
}

int cvOneDMaterialManager::AddNewMaterialLinear(double density, double dynamicViscosity,
                                                double profile_exponent, double pRef,
                                                double EHR){
  cvOneDMaterialLinear* linearmat = new cvOneDMaterialLinear();
  linearmat->SetDensity(density);
  linearmat->SetDynamicViscosity(dynamicViscosity);
  linearmat->SetProfileExponent(profile_exponent);
  linearmat->SetReferencePressure(pRef);
  linearmat->SetEHR(EHR,pRef);
  return cvOneDGlobal::gMaterialManager->AddNewMaterial(MaterialType_MATERIAL_LINEAR,(cvOneDMaterial*)linearmat);
}

// caller must deallocate material instance to avoid memory leak
cvOneDMaterial* cvOneDMaterialManager::GetNewInstance(int matID){
  if (types[matID] == MaterialType_MATERIAL_OLUFSEN) {
    cvOneDMaterialOlufsen* olfmat = new cvOneDMaterialOlufsen();
  //  printf("In GetNewInstance cvOneDMaterialOlufsen is called  matID=%i \n",matID);
    *olfmat = *((cvOneDMaterialOlufsen*)(materials[matID]));
  //  printf("In GetNewInstance cvOneDMaterialOlufsen* materials is called \n");
    return (cvOneDMaterial*)olfmat;
  }else if (types[matID] == MaterialType_MATERIAL_LINEAR) {
    cvOneDMaterialLinear* linearmat = new cvOneDMaterialLinear();
    *linearmat = *((cvOneDMaterialLinear*)(materials[matID]));
    return (cvOneDMaterial*)linearmat;
  }else{
    return NULL;
  }
}
