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
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include "vtkExecutive.h"
#include <sstream>
#include <vtkPointSet.h>
#include "hdf5.h"


#include <QApplication>
#include <QStyle>

#include "pqApplicationCore.h"
#include "pqObjectBuilder.h"
#include "pqServer.h"
#include "pqServerManagerModel.h"
#include "pqUndoStack.h"

class pvESSI : public vtkUnstructuredGridAlgorithm
{
public:
  vtkTypeMacro(pvESSI,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  void Plot_Node_Mesh(int x){ if(x) Display_Node_Mesh=1;	else Display_Node_Mesh =0; }
  void Plot_Gauss_Mesh(int x){ if(x){/**this->Get_Gauss_Mesh(UGrid_Gauss_Mesh);**/ Display_Gauss_Mesh=1;}  else Display_Gauss_Mesh =0; }
  void Plot_All_Mesh(int x){ if(x) Display_All_Mesh=1;  else Display_All_Mesh =0; }
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
  int ProcessRequest(vtkInformation *, vtkInformationVector ** , vtkInformationVector * );
  static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);
  
  //////////////////////// Important Variables /////////////////////////////////////////

  int TimeStep, Number_Of_Time_Steps;;
  int Node_Mesh_Current_Time, Gauss_Mesh_Current_Time;
  int Node_Mesh_Previous_Time, Gauss_Mesh_Previous_Time; 

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
  std::map<int,std::vector<int> > ESSI_to_VTK_Connectivity;
  std::map<double,int > Time_Map;
  pvESSI(const pvESSI&);  // Not implemented.
  void operator=(const pvESSI&);  // Not implemented.
  void set_VTK_To_ESSI_Elements_Connectivity();
  void initialize();
  void Build_Maps();
  void Build_Time_Map();
  void Build_Gauss_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Gauss_Mesh);
  void Build_Node_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh);
  void Build_All_Attributes(vtkSmartPointer<vtkUnstructuredGrid> All_Mesh); // need to implement
  void Build_Delaunay3D_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Mesh);
  void Build_ProbeFilter_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Probe_Input, int probe_type);  // Probing variables at gauss nodes from node mesh
  double *Time; 
 
  char* FileName;

  /******************************************* Mesh ******************************************/
  vtkSmartPointer<vtkUnstructuredGrid> UGrid_Node_Mesh;          // Mesh with nodes
  vtkSmartPointer<vtkUnstructuredGrid> UGrid_Gauss_Mesh;         // Mesh with only gauss points
  vtkSmartPointer<vtkUnstructuredGrid> UGrid_All_Mesh;           // Mesh with both gauss points and nodes  // Not implemented
  vtkSmartPointer<vtkUnstructuredGrid> UGrid_Current_Node_Mesh;  // Contains mesh with data attributes 
  vtkSmartPointer<vtkUnstructuredGrid> UGrid_Current_Gauss_Mesh; // Contains mesh with data attributes   
  vtkSmartPointer<vtkUnstructuredGrid> UGrid_Current_All_Mesh;   // Contains mesh with data attributes     // Not implemented

  int Number_of_Elements, Number_of_Nodes, Number_of_Gauss_Nodes;
  int Pseudo_Number_of_Elements, Pseudo_Number_of_Nodes;
  int Display_Node_Mesh,Display_Gauss_Mesh,Display_All_Mesh,Whether_Node_Mesh_Build, Whether_Gauss_Mesh_Build, Whether_All_Mesh_Build;


  void Get_Node_Mesh(vtkSmartPointer<vtkUnstructuredGrid> UGrid_Node_Mesh); 	  // Building the node mesh skeleton
  void Get_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> UGrid_Gauss_Mesh); 	// Building the gauss mesh skeleton
  void Get_All_Mesh(vtkSmartPointer<vtkUnstructuredGrid> All_Mesh);             // Building the node as well as gauss mesh skeleton // still to implement
  void SetMetaArrays( vtkSmartPointer<vtkFloatArray> &vtk_Generalized_Displacements, vtkSmartPointer<vtkFloatArray> &vtk_Generalized_Velocity, 
    vtkSmartPointer<vtkFloatArray> &vtk_Generalized_Acceleration, vtkSmartPointer<vtkFloatArray> &Elastic_Strain_Tensor, vtkSmartPointer<vtkFloatArray> &Plastic_Strain_Tensor, 
    vtkSmartPointer<vtkFloatArray> &Stress_Tensor, vtkSmartPointer<vtkFloatArray> &Material_Properties, vtkSmartPointer<vtkFloatArray> &Total_Energy, 
    vtkSmartPointer<vtkFloatArray> &Incremental_Energy, vtkSmartPointer<vtkIntArray> &Node_Label, vtkSmartPointer<vtkIntArray> &Element_Label );

  /*******************************************************************************************/

};
 
#endif
