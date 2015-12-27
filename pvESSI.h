#ifndef __pvESSI_h
#define __pvESSI_h
 
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkQuadratureSchemeDictionaryGenerator.h"
#include "vtkInformationQuadratureSchemeDefinitionVectorKey.h"
#include "vtkInformationStringKey.h"
#include "vtkFloatArray.h"
#include "vtkSmartPointer.h"
#include <map>
#include <vector>
#include <iostream>
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkQuadratureSchemeDefinition.h"
#include "vtkQuadraturePointInterpolator.h"
#include "vtkQuadraturePointsGenerator.h"
#include "vtkQuadratureSchemeDictionaryGenerator.h"
#include "vtkInformationQuadratureSchemeDefinitionVectorKey.h"
#include "vtkIdTypeArray.h"
#include "vtkInformationVector.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkLookupTable.h"
#include "H5Cpp.h"
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include "vtkExecutive.h"
#include <sstream>

class pvESSI : public vtkUnstructuredGridAlgorithm
{
public:
  vtkTypeMacro(pvESSI,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  void PlotGaussMesh(int x){
  	if(x) {
  		std::cout << "Displaying Gauss Mesh" <<std::endl;
  		Display_Gauss_Mesh=1;
  	}
  	else{
  		std::cout << "Displaying General Mesh" <<std::endl;
  		Display_Gauss_Mesh =0;
  	}
  }

  void PlotGeneralMesh(int x){
  	if(x) {
  		std::cout << "Displaying General Mesh" <<std::endl;
  		Display_General_Mesh =1;
  	}
  	else{
  		std::cout << "Displaying Gauss Mesh" <<std::endl;
  		Display_General_Mesh =0;
  	}
  }

  void PrintX(int x){}

  // static vtkInformationQuadratureSchemeDefinitionVectorKey* DICTIONARY();
  // static vtkInformationStringKey* QUADRATURE_OFFSET_ARRAY_NAME();

  // void GetNumberOfCellArrays();
 
  static pvESSI *New();
 
  // Description:
  // Specify file name of the .abc file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
  vtkSetMacro(TimeStep, int);
  vtkGetMacro(TimeStep, int);
  int GetNumberOfTimeSteps(){
    return this->Number_Of_Time_Steps;
  }
 
protected:
  pvESSI();
  ~pvESSI(){}
 
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestInformation( vtkInformation *, vtkInformationVector **, vtkInformationVector* );
  static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);
  
  //////////////////////// Important Variables /////////////////////////////////////////

  int TimeStep;
  int Current_Time; 
  int Number_Of_Time_Steps;
  // vtkTimeStamp ReadMTime;
  // int ReadError;
  int Time_Step_Range[2];

// // incrementally fine-tuned progress updates.
//   virtual void GetProgressRange(float* range);
//   virtual void SetProgressRange(const float range[2], int curStep, int numSteps);
//   virtual void SetProgressRange(const float range[2], int curStep, const float* fractions);
//   virtual void UpdateProgressDiscrete(float progress);
//   float ProgressRange[2];

  vtkDataObject* GetCurrentOutput();
  vtkInformation* GetCurrentOutputInformation();
 
private:

  std::map<int,int> ESSI_to_VTK_Element;
  std::map< int,std::vector<int> > ESSI_to_VTK_Connectivity;
  pvESSI(const pvESSI&);  // Not implemented.
  void operator=(const pvESSI&);  // Not implemented.
  void set_VTK_To_ESSI_Elements_Connectivity();
  void Build_Gauss_Attributes();
  void Build_Node_Attributes();
  int Energy_Database_Status;
  float *Prev_Total_Energy_Database;
 
  char* FileName;

  /******************************************* Mesh ******************************************/
  vtkSmartPointer<vtkUnstructuredGrid> UGrid_Mesh;
  vtkSmartPointer<vtkUnstructuredGrid> UGrid_Gauss_Mesh;
  vtkSmartPointer<vtkUnstructuredGrid> *UGrid_All_Mesh;
  int Number_of_Elements, Number_of_Nodes, Number_of_Gauss_Nodes;
  int Display_Gauss_Mesh=0, Display_General_Mesh =1, whetehr_general_mesh_build=0, whetehr_gauss_mesh_build=0;

  void GetGeneralMesh(); 	// Building the general mesh skeleton
  void GetGaussMesh(); 		// Building the gauss mesh skeleton
  void SetMetaArrays( vtkSmartPointer<vtkFloatArray> &vtk_Generalized_Displacements, vtkSmartPointer<vtkFloatArray> &vtk_Generalized_Velocity, vtkSmartPointer<vtkFloatArray> &vtk_Generalized_Acceleration, vtkSmartPointer<vtkFloatArray> &Elastic_Strain_Tensor, vtkSmartPointer<vtkFloatArray> &Plastic_Strain_Tensor, vtkSmartPointer<vtkFloatArray> &Stress_Tensor, vtkSmartPointer<vtkFloatArray> &Material_Properties, vtkSmartPointer<vtkFloatArray> &Total_Energy, vtkSmartPointer<vtkFloatArray> &Incremental_Energy);


  /*******************************************************************************************/

};
 
#endif
