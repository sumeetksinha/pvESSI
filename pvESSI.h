#ifndef __pvESSI_h
#define __pvESSI_h
 
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkQuadratureSchemeDictionaryGenerator.h"
#include "vtkInformationQuadratureSchemeDefinitionVectorKey.h"
#include "vtkInformationIntegerKey.h"
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
#include <boost/regex.hpp>
#include <vtkAppendFilter.h>
#include <QApplication>
#include <QStyle>
#include "vtkDataArraySelection.h"

#include "pqApplicationCore.h"
#include "pqObjectBuilder.h"
#include "pqServer.h"
#include "pqServerManagerModel.h"
#include "pqUndoStack.h"

class pvESSI : public vtkUnstructuredGridAlgorithm{

public:
  vtkTypeMacro(pvESSI,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  void Enable_Gauss_To_Node_Interpolation(int x){ if(!(x and Enable_Gauss_To_Node_Interpolation_Flag))this->Modified();  Enable_Gauss_To_Node_Interpolation_Flag=false; if(x) Enable_Gauss_To_Node_Interpolation_Flag=true; }
  void Show_Gauss_Mesh(int x){ if(!(x and Show_Gauss_Mesh_Flag))this->Modified(); Show_Gauss_Mesh_Flag=false; if(x) Show_Gauss_Mesh_Flag=true;}
  void Build_pvESSI_Folder(int x){ if(!(x and Enable_Building_of_Maps_Flag))this->Modified(); Enable_Building_of_Maps_Flag=true; if(x) Enable_Building_of_Maps_Flag=false; Enable_Initialization_Flag=true;}
  void Enable_Relative_Displacement(int x ){if(!(x and Enable_Relative_Displacement_Flag))this->Modified();  Enable_Relative_Displacement_Flag=false; if(x) Enable_Relative_Displacement_Flag=true;}
  void Enable_Displacement_Probing(int x ){ if(!(x and Enable_Displacement_Probing_Flag))this->Modified();   Enable_Displacement_Probing_Flag=false; if(x) Enable_Displacement_Probing_Flag=true;}
  void Reference_Displacement_Index(int x ){if(x!=Reference_Displacement_Index_Flag)this->Modified();   Reference_Displacement_Index_Flag=x;}
  void Enable_Physical_Node_Group_Selection(int x ){ if(!(x and Enable_Physical_Node_Group_Selection_Flag))this->Modified();  Enable_Physical_Node_Group_Selection_Flag=false; if(x) Enable_Physical_Node_Group_Selection_Flag=true;}
  void Enable_Physical_Element_Group_Selection(int x ){ if(!(x and Enable_Physical_Element_Group_Selection_Flag))this->Modified();  Enable_Physical_Element_Group_Selection_Flag=false; if(x) Enable_Physical_Element_Group_Selection_Flag=true;}

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
    return this->Number_of_Time_Steps;
  }

  /*************** Physical Element Groups ****************************/
  void SetPhysicalElementGroupArrayStatus(const char* name, int status);
  const char* GetPhysicalElementGroupArrayName(int index);
  int GetNumberOfPhysicalElementGroupArrays();
  int GetPhysicalElementGroupArrayStatus(const char* name);

  /*************** Physical Node Groups ****************************/
  void SetPhysicalNodeGroupArrayStatus(const char* name, int status);
  const char* GetPhysicalNodeGroupArrayName(int index);
  int GetNumberOfPhysicalNodeGroupArrays();
  int GetPhysicalNodeGroupArrayStatus(const char* name);

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
  int Number_of_Time_Steps;           // Stores the Number of TimeSteps in the analysis
  int Number_of_Sub_Steps;            // Stores the Number of Substeps in the analysis
  int Node_Mesh_Current_Time;         // Stores the currengt time step of Node Mesh
  int Gauss_Mesh_Current_Time;        // Stores the currengt time step of Gauss Mesh

  /************************************ Element Info Parameters *********************************/
  int Pseudo_Number_of_Elements;      // Stores the max node tag in each domain
  int Number_of_Elements;             // Stores the Number of Elements in each domain
  int Number_of_Gauss_Points;         // Stores the Number of gauss points in each domain
  int Number_of_Connectivity_Nodes;   // Stores the Number of connectivity nodes in reach domain

  /********************************** Node Info parameters **************************************/
  int Number_of_Nodes;                // Stores the numbr of nodes in each domain
  int Pseudo_Number_of_Nodes;         // Stores the max Node Tag in each demain 
  int Number_of_Constrained_Dofs;     // Number of constrained dofs

  /********************************* Output Class_DESC_ENCODING Format *************************/
  int ELE_TAG_DESC_ENCODING;

  /********************************** Model Info *************************************************/
  int Number_of_Processes_Used;       // Number of Processes used
  int Process_Number;
  bool single_file_visualization_mode; 
  
  /*************************** Visualization Parameters *****************************************/
  int Build_Map_Status;               // Whether Map is Build
  bool *Whether_Node_Mesh_build;      // Whether Node Mesh Build For Domains
  bool *Whether_Gauss_Mesh_build;     // Whether Node Mesh Build For Domains
  bool Enable_Initialization_Flag;    // Whether Initialization Done
  int piece_no;                       // Piece no or Processor no
  int num_of_pieces;                  // total number of pieces or processors
  int Number_of_Strain_Strain_Info;   // Total Number_of_Stress_Strain_Data_at_Nodes
  int domain_no                   ;   // domain no of the mesh
  bool enable_support_reactions;
  bool eigen_mode_on;                 // enable eigen analysis mode

  ///////////////////////////// HDF5 ID /////////////////////////////////////////////////////// 

  /***************** File_id **********************************/
  hid_t id_File;

  /***************** Time Steps *******************************/
  hid_t id_time;
  hid_t id_Number_of_Time_Steps;

  /***************** Model Info *******************************/
  hid_t id_Model_group; 
  hid_t id_Number_of_Elements;
  hid_t id_Number_of_Nodes;
  hid_t id_Number_of_Gauss_Points;
  hid_t id_Number_of_Processes_Used;
  hid_t id_Process_Number;

  /**************** Element Info ******************************/
  hid_t id_Elements_group;
  hid_t id_Class_Tags;
  hid_t id_Element_Class_Desc;
  hid_t id_Connectivity;
  hid_t id_Gauss_Point_Coordinates;
  hid_t id_Index_to_Connectivity;
  hid_t id_Material_Tags;
  hid_t id_Element_Outputs;
  hid_t id_Gauss_Outputs;

  /**************** Node Info ******************************/
  hid_t id_Nodes_group;
  hid_t id_Constrained_DOFs;
  hid_t id_Constrained_Nodes;
  hid_t id_Coordinates;
  hid_t id_Generalized_Displacements;
  hid_t id_Support_Reactions;
  hid_t id_Number_of_DOFs;

  /**************** Maps ***********************************/
  hid_t id_pvESSI; 
  hid_t id_Element_Map;
  hid_t id_Node_Map;
  hid_t id_Inverse_Node_Map;
  hid_t id_Inverse_Element_Map;
  hid_t id_Number_of_Elements_Shared;
  hid_t id_Number_of_Gauss_Elements_Shared;

  /*************** Partition Info *************************/
  hid_t id_Node_Partition;
  hid_t id_Element_Partition;

  /*************** Field at Nodes ***************************/
  hid_t id_Field_at_Nodes_group;
  hid_t id_Stress_and_Strain;
  hid_t id_Whether_Stress_Strain_Build;
  hid_t id_Energy;                            // Not implemented
  hid_t id_Whether_Energy_Build;              // Not implemented

  /************** Substep Outputs ***************************/
  hid_t id_Number_of_Iterations;
  hid_t id_Iterative_Generalized_Displacements;
  hid_t id_Iterative_Element_Outputs;
  hid_t id_Iterative_Gauss_Outputs;

  /************* Eigen Mode Analysis ***********************/
  hid_t id_Eigen_Mode_Analysis;
  hid_t id_frequencies;
  hid_t id_modes;
  hid_t id_number_of_modes;
  hid_t id_periods;
  hid_t id_values;

  /************** Physical Groups ***************************/
  hid_t id_Physical_Element_Groups;
  hid_t id_Physical_Node_Groups;

  /************** General Variable **************************/
  hid_t DataSpace;
  hid_t DataSet;
  hid_t Group; 
  hid_t MemSpace;

  /************** Visualization Control Variables **********/
  bool Enable_Gauss_To_Node_Interpolation_Flag;
  bool Enable_Building_of_Maps_Flag;
  bool Enable_Relative_Displacement_Flag;
  int Reference_Displacement_Index_Flag;
  bool Show_Gauss_Mesh_Flag;
  bool Enable_Displacement_Probing_Flag;
  bool Whether_Physical_Group_Info_build;
  bool Enable_Physical_Element_Group_Selection_Flag;
  bool Enable_Physical_Node_Group_Selection_Flag;
  vtkSmartPointer<vtkDataArraySelection> Physical_Node_Group;
  vtkSmartPointer<vtkDataArraySelection> Physical_Element_Group;



  hsize_t  dims1_out[1], dims2_out[2];
  hsize_t  dims3[3],     dims2[2];
  hsize_t  maxdims3[3],  maxdims2[2];

  hsize_t  block3[3],  block2[2],  block1[1];
  hsize_t  count3[3],  count2[2],  count1[1];
  hsize_t offset3[3], offset2[2], offset1[1];
  hsize_t stride3[3], stride2[2], stride1[1];

  herr_t status;

  hsize_t index_i,index_j,index_k;

  int node_no, element_no;

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

static herr_t op_func (hid_t loc_id, const char *name, const H5O_info_t *info, void *operator_data);


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
  void Domain_Initializer(int Piece_No); 
  void Close_File();
  void Build_Maps();
  void Build_Meta_Array_Map();
  void Build_Gauss_To_Node_Interpolation_Map();
  void Build_Gauss_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Gauss_Mesh, int Node_Mesh_Current_Time);
  void Build_Node_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh, int Gauss_Mesh_Current_Time);
  void Build_Eigen_Modes_Node_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh, int Current_Time);
  void Build_Delaunay3D_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Mesh);
  void Build_ProbeFilter_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Probe_Input, int Current_Time); // Probing variables at gauss nodes from node mesh
  void Build_Stress_Field_At_Nodes_v2(vtkSmartPointer<vtkUnstructuredGrid> Gauss_Mesh, int Node_Mesh_Current_Time);
  void Build_Stress_Field_At_Nodes(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh, int Node_Mesh_Current_Time);
  void Build_Physical_Element_Group_Mesh(vtkSmartPointer<vtkUnstructuredGrid> NodeMesh);
  void Build_VTK_Element_Selection_mesh();
  std::string GetSourceFile(std::string filename);


  double *Time;            // Holds the time vector for the domain
  char* FileName;          // Holds the filename
  int* Element_Desc_Array; // Holds the element decription array 

  /******************************************* Mesh ******************************************/  
  vtkSmartPointer<vtkUnstructuredGrid> *UGrid_Node_Mesh;                 // Contains the mesh of all domains
  vtkSmartPointer<vtkUnstructuredGrid> *UGrid_Gauss_Mesh;                // Contains the mesh of all domains
  vtkSmartPointer<vtkUnstructuredGrid> *UGrid_Current_Node_Mesh;         // Contains mesh with data attributes 
  vtkSmartPointer<vtkUnstructuredGrid> *UGrid_Current_Gauss_Mesh;        // Contains mesh with data attributes 

  /****************************************** Physical Groups*********************************/
  static std::vector<std::string> Physical_Group_Container;
  /******************************* Meta Data Arrays ********************************************/

  // Meta Data Arrays
  vtkSmartPointer<vtkFloatArray> Generalized_Displacements;
  vtkSmartPointer<vtkIntArray>   Boundary_Conditions;
  vtkSmartPointer<vtkFloatArray> Support_Reactions;
  vtkSmartPointer<vtkFloatArray> Generalized_Forces;
  vtkSmartPointer<vtkFloatArray> Generalized_Velocity;
  vtkSmartPointer<vtkFloatArray> Generalized_Acceleration;

  // Stress-Strain 
  vtkSmartPointer<vtkFloatArray> Elastic_Strain;
  vtkSmartPointer<vtkFloatArray> Plastic_Strain;
  vtkSmartPointer<vtkFloatArray> Stress;

  // Stress-Strain Invariants
  vtkSmartPointer<vtkFloatArray> q;
  vtkSmartPointer<vtkFloatArray> p;
  vtkSmartPointer<vtkFloatArray> Plastic_Strain_q;
  vtkSmartPointer<vtkFloatArray> Plastic_Strain_p;

  // Tags
  vtkSmartPointer<vtkIntArray>   Material_Tag;
  vtkSmartPointer<vtkIntArray>   Node_Tag;
  vtkSmartPointer<vtkIntArray>   Element_Tag;
  vtkSmartPointer<vtkIntArray>   Class_Tag;

  // Energy 
  vtkSmartPointer<vtkIntArray>   Total_Energy;
  vtkSmartPointer<vtkIntArray>   Incremental_Energy;  

  // Partition Information 
  vtkSmartPointer<vtkIntArray>   Partition_Info; // For Elements


  /***************************** Mesh Building Functions ******************************************/

  void Get_Node_Mesh(vtkSmartPointer<vtkUnstructuredGrid> UGrid_Node_Mesh); 	    // Building the node mesh skeleton
  void Get_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> UGrid_Gauss_Mesh); 	  // Building the gauss mesh skeleton
  void Merge_Mesh(int start, int end, vtkSmartPointer<vtkUnstructuredGrid> Mesh); // Merge domain mesh  when necessary
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
