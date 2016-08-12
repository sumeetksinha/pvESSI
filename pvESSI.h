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

class pvESSI : public vtkUnstructuredGridAlgorithm{

public:
  vtkTypeMacro(pvESSI,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  void Plot_Node_Mesh(int x){ /*if(x) cout >>;	else Display_Node_Mesh =0;*/ }
  void Plot_Gauss_Mesh(int x){ /*if(x){/**this->Get_Gauss_Mesh(UGrid_Gauss_Mesh);**//* Display_Gauss_Mesh=1;}  else Display_Gauss_Mesh =0;*/ }
  void Plot_All_Mesh(int x){/* if(x) Display_All_Mesh=1;  else Display_All_Mesh =0; */}
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

  /************************************** Time parameters ***************************************/
  int TimeStep;                 
  int Time_Step_Range[2];             // Range of Time Steps
  int Number_Of_Time_Steps;           // Stores the Number of TimeSteps in the analysis
  int Node_Mesh_Current_Time;         // Stores the currengt time step of Node Mesh
  int Gauss_Mesh_Current_Time;        // Stores the currengt time step of Gauss Mesh

  /************************************ Element Info Parameters *********************************/
  int Number_of_Gauss_Nodes;          // Stores the Number of Gauss Nodes in Model
  int Pseudo_Number_of_Elements;      // Stores the Number of Elements in ESSI ( max_ESSI_Element_Tag )
  int Number_of_Elements;             // Stores the Number of Elements in the model

  /********************************** Node Info parameters **************************************/
  int Number_of_Nodes;                // Stores the Number of Nodes in Model
  int Pseudo_Number_of_Nodes;         // Stores the Number of Nodes in ESSI ( max_ESSI_Node_Tag ) 

  /********************************** Model Info *************************************************/
  int Number_of_Processes_Used;       // Number of Processes used
  int Process_Number;
  
  /*************************** Visualization Parameters *****************************************/
  int Display_Node_Mesh;              // Whether One Wants to display Node Mesh
  int Display_Gauss_Mesh;             // Whether one wants to display gauss mesh
  int Whether_Node_Mesh_Build;        // Whether node_mesh_build
  int Whether_Gauss_Mesh_Build;       // Whether gauss mesh build
  int Build_Map_Status;               // Whether Map is Build
  bool Enable_Gauss_Mesh=false;       // Enable Gauss Mesh  
  int EXTENT[6];                      // Extent in int
  float Model_Bounds[6];              // Model Bound in float

  ///////////////////////////// HDF5 ID /////////////////////////////////////////////////////// 

  /***************** File_id **********************************/
  hid_t id_File;

  /***************** Time Steps *******************************/
  hid_t id_time;
  hid_t id_Number_of_Time_Steps;

  /***************** Model Info *******************************/
  hid_t id_Model_group; 
  hid_t id_Model_Bounds;
  hid_t id_Whether_Maps_Build; 
  hid_t id_Number_of_Elements;
  hid_t id_Number_of_Nodes;
  hid_t id_Number_of_Processes_Used;
  hid_t id_Process_Number;

  /**************** Element Info ******************************/
  hid_t id_Elements_group;
  hid_t id_Class_Tags;
  hid_t id_Connectivity;
  hid_t id_Element_types;
  hid_t id_Gauss_Point_Coordinates;
  hid_t id_Index_to_Connectivity;
  hid_t id_Index_to_Gauss_Point_Coordinates;
  hid_t id_Index_to_Outputs;
  hid_t id_Material_tags;
  hid_t id_Number_of_Gauss_Points;
  hid_t id_Element_Number_of_Nodes;
  hid_t id_Number_of_Output_Fields;
  hid_t id_Outputs;
  hid_t id_Substep_Outputs;

  /**************** Node Info ******************************/
  hid_t id_Nodes_group;
  hid_t id_Constrained_DOFs;
  hid_t id_Constarined_Nodes;
  hid_t id_Coordinates;
  hid_t id_Generalized_Displacements;
  hid_t id_Index_to_Coordinates;
  hid_t id_Index_to_Generalized_Displacements;
  hid_t id_Number_of_DOFs;

  /**************** Maps ***********************************/
  hid_t id_Maps_group;
  hid_t id_Element_Map;
  hid_t id_Node_Map;
  hid_t id_Inverse_Node_Map;
  hid_t id_Inverse_Element_Map;
  hid_t id_Number_of_Elements_Shared;
  hid_t id_Number_of_Gauss_Elements_Shared;

  /*************** Field at Nodes ***************************/
  hid_t id_Field_at_Nodes_group;
  hid_t id_Stress_and_Strain;
  hid_t id_Whether_Stress_Strain_Build;
  hid_t id_Energy;                            // Not implemented
  hid_t id_Whether_Energy_Build;              // Not implemented


  /************** General Variable **************************/
  hid_t DataSpace;
  hid_t DataSet;
  hid_t Group; 
  hid_t MemSpace;

  hsize_t  dims1_out[1], dims2_out[2];
  hsize_t  dims3[3],     dims2[2];
  hsize_t  maxdims3[3],  maxdims2[2];

  hsize_t  block3[3],  block2[2],  block1[1];
  hsize_t  count3[3],  count2[2],  count1[1];
  hsize_t offset3[3], offset2[2], offset1[1];
  hsize_t stride3[3], stride2[2], stride1[1];

  herr_t status;

  hsize_t index_i,index_j,index_k;

  int Int_Variable_1, Int_Variable_2, Int_Variable_3, Int_Variable_4, index;
  int node_no, element_no;
  float Float_Variable_1, Float_Variable_2;

  /************* Hdf5 function ******************************/
  void HDF5_Read_INT_Array_Data(hid_t id_DataSet,
                               int rank,
                               hsize_t *data_dims,
                               hsize_t *offset,
                               hsize_t *stride,
                               hsize_t *count,
                               hsize_t *block,
                               int* data);

  void HDF5_Write_INT_Array_Data(hid_t id_DataSet,
                                 int rank,
                                 hsize_t *data_dims,
                                 hsize_t *offset,
                                 hsize_t *stride,
                                 hsize_t *count,
                                 hsize_t *block,
                                 int* data);

  void HDF5_Read_FLOAT_Array_Data(hid_t id_DataSet,
                                   int rank,
                                   hsize_t *data_dims,
                                   hsize_t *offset,
                                   hsize_t *stride,
                                   hsize_t *count,
                                   hsize_t *block,
                                   float* data);

  void HDF5_Read_DOUBLE_Array_Data(hid_t id_DataSet,
                                   int rank,
                                   hsize_t *data_dims,
                                   hsize_t *offset,
                                   hsize_t *stride,
                                   hsize_t *count,
                                   hsize_t *block,
                                   double* data);

  void HDF5_Write_FLOAT_Array_Data(hid_t id_DataSet,
                                   int rank,
                                   hsize_t *data_dims,
                                   hsize_t *offset,
                                   hsize_t *stride,
                                   hsize_t *count,
                                   hsize_t *block,
                                   float* data);

  void HDF5_Write_DOUBLE_Array_Data(hid_t id_DataSet,
                                   int rank,
                                   hsize_t *data_dims,
                                   hsize_t *offset,
                                   hsize_t *stride,
                                   hsize_t *count,
                                   hsize_t *block,
                                   double* data);

  void HDF5_Read_STRING_Array_Data(hid_t id_DataSet,
                                   int rank,
                                   hsize_t *data_dims,
                                   hsize_t *offset,
                                   hsize_t *stride,
                                   hsize_t *count,
                                   hsize_t *block,
                                   int* data);

  void Append_Number_of_Elements_Shared(int message);

  /************************* Some common Variables ********************************/

  // vtkTimeStamp ReadMTime;
  // int ReadError;

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
  std::map<std::string,int > Meta_Array_Map;
  std::map<int,double**> Gauss_To_Node_Interpolation_Map;

  pvESSI(const pvESSI&);  // Not implemented.
  void operator=(const pvESSI&);  // Not implemented.
  void set_VTK_To_ESSI_Elements_Connectivity();
  void Initialize();
  void Build_Time_Map();
  void Step_Initializer(int Piece_No); 
  void Close_File();
  void Build_Maps();
  void Build_Meta_Array_Map();
  void Build_Gauss_To_Node_Interpolation_Map();
  void Build_Gauss_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Gauss_Mesh, int Node_Mesh_Current_Time);
  void Build_Node_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh, int Gauss_Mesh_Current_Time);
  void Build_Delaunay3D_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Mesh);
  void Build_ProbeFilter_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Probe_Input, int probe_type);  // Probing variables at gauss nodes from node mesh
  void Build_Stress_Field_At_Nodes_v2(vtkSmartPointer<vtkUnstructuredGrid> Gauss_Mesh, int Node_Mesh_Current_Time);
  void Build_Stress_Field_At_Nodes(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh, int Node_Mesh_Current_Time);
  double *Time; 
 
  char* FileName;

  /******************************************* Mesh ******************************************/
  vtkSmartPointer<vtkUnstructuredGrid> UGrid_Node_Mesh;          // Mesh with nodes
  vtkSmartPointer<vtkUnstructuredGrid> UGrid_Gauss_Mesh;         // Mesh with only gauss points
  vtkSmartPointer<vtkUnstructuredGrid> UGrid_Current_Node_Mesh;  // Contains mesh with data attributes 
  vtkSmartPointer<vtkUnstructuredGrid> UGrid_Current_Gauss_Mesh; // Contains mesh with data attributes   
  // vtkSmartPointer<vtkUnstructure

  ///////////////////////////// HDF5 ID /////////////////////////////////////////////////////// dGrid> UGrid_Current_All_Mesh;   // Contains mesh with data attributes     // Not implemented

  // Meta Da
  vtkSmartPointer<vtkFloatArray> Generalized_Displacements;
  vtkSmartPointer<vtkFloatArray> Generalized_Velocity;
  vtkSmartPointer<vtkFloatArray> Generalized_Acceleration;
  vtkSmartPointer<vtkFloatArray> Elastic_Strain;
  vtkSmartPointer<vtkFloatArray> Plastic_Strain;
  vtkSmartPointer<vtkFloatArray> Stress;
  vtkSmartPointer<vtkIntArray> Material_Tag;
  vtkSmartPointer<vtkIntArray> Total_Energy;
  vtkSmartPointer<vtkIntArray> Incremental_Energy;
  vtkSmartPointer<vtkIntArray> Node_Tag;
  vtkSmartPointer<vtkIntArray> Element_Tag;


  void Get_Node_Mesh(vtkSmartPointer<vtkUnstructuredGrid> UGrid_Node_Mesh); 	  // Building the node mesh skeleton
  void Get_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> UGrid_Gauss_Mesh); 	// Building the gauss mesh skeleton
  void Set_Meta_Array(int Meta_Data_Id );

  void Build_Inverse_Matrices();
  void Build_Brick_Coordinates();
  double Brick_8_Gauss_Coordinates[8][3];
  double Brick_20_Gauss_Coordinates[20][3];
  double Brick_27_Gauss_Coordinates[27][3];
  double Brick_Coordinates[27][3];
  double **Twenty_Node_Brick_Inverse;
  double **Twenty_Seven_Node_Brick_Inverse;
  double **Eight_Node_Brick_Inverse;
  /*******************************************************************************************/

};
 
#endif
