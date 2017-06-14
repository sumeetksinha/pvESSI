#include "pvESSI.h"
#include <stdlib.h>
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkQuadratureSchemeDefinition.h"
#include "vtkQuadraturePointInterpolator.h"
#include "vtkQuadraturePointsGenerator.h"
#include "vtkQuadratureSchemeDictionaryGenerator.h"
#include "vtkInformationQuadratureSchemeDefinitionVectorKey.h"
#include "vtkIdTypeArray.h"
#include "vtkInformationVector.h"
#include "vtkDataArray.h"
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
#include "vtkDelaunay3D.h"
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include "vtkExecutive.h"
#include "vtkPVView.h"
#include <sstream>
#include "hdf5.h"
#include "vtkProbeFilter.h"
#include <vtkDelaunay3D.h>
#include "vtkSelectionNode.h"
#include "vtkExtractSelection.h"
#include "vtkSelection.h"
#include "vtkGaussianSplatter.h"


// #define DEBUG_MODE 1
#define PRINT_FUNCTION(function)  	cout <<"\n====================================\n"<<function<<"\n====================================\n\n";

// LINK_LIBRARIES(hdf5_cpp hdf5 )

//============================================================================
// Formula to get Number of Nodes, Number of Gauss, Number of Element Output
// for a given class_tag_desc and given class encoding 
//============================================================================/
#define NUMBER_OF_NODES(class_tag_desc)             (class_tag_desc/(int)pow(10.0,(9-(ELE_TAG_DESC_ENCODING/(int)pow(10.0,7))%100%10)))%((int)pow(10.0,((ELE_TAG_DESC_ENCODING/(int)pow(10.0,7))%100%10   - (ELE_TAG_DESC_ENCODING/(int)pow(10.0,8))%100%10  + 1)))
#define NUMBER_OF_GAUSS(class_tag_desc)             (class_tag_desc/(int)pow(10.0,(9-(ELE_TAG_DESC_ENCODING/(int)pow(10.0,5))%100%10)))%((int)pow(10.0,((ELE_TAG_DESC_ENCODING/(int)pow(10.0,5))%100%10   - (ELE_TAG_DESC_ENCODING/(int)pow(10.0,6))%100%10  + 1))) 
#define NUMBER_OF_ELEMENT_OUTPUTS(class_tag_desc)   (class_tag_desc/(int)pow(10.0,(9-(ELE_TAG_DESC_ENCODING/(int)pow(10.0,3))%100%10)))%((int)pow(10.0,((ELE_TAG_DESC_ENCODING/(int)pow(10.0,3))%100%10   - (ELE_TAG_DESC_ENCODING/(int)pow(10.0,4))%100%10  + 1))) 
#define ELEMENT_CATEGORY(class_tag_desc)   (class_tag_desc/(int)pow(10.0,8))


std::vector<std::string> pvESSI::Physical_Group_Container;
vtkStandardNewMacro(pvESSI);

//============================================================================
// pvESSI  constructor called when the plugin is loaded in paraview 
//      
// * All whether_dataset_arrays musst be initialized to false. 
// * Whether physical group information is build is set to false.
// * Iniitialization flag is also set to false as it needs to be initialized.
//   for the first time.
// * Set the number of output ports as 1 and input ports as 0.
// * Build VTK_to_ESSI_Connectivity and Meta_Array Maps.
// * Build Inverse_Matrices and Gauss_To_Node_Interpolation_Map.
//============================================================================
pvESSI::pvESSI(){ 

	this->FileName = NULL;                                   // set file to NULL 
	this->eigen_mode_on = false;                             // set eigen mode to false
	this->Whether_Node_Mesh_Array_Initialized = false;       // set whether_node_mesh_array initialized to be false 
	this->Whether_Node_Mesh_Attributes_Initialized = false;  // set whether_node_mesh_attribute_array initialized to be false 
	this->Whether_Gauss_Mesh_Array_Initialized = false;      // set whether_gauss_mesh_array initialized to be false 
	this->Whether_Gauss_Mesh_Attributes_Initialized = false; // set whether_gauss_mesh_attribute_array initialized to be false 
	this->Whether_Node_Mesh_Stress_Attributes_Initialized=false;

	this->Enable_Initialization_Flag=true;              	 // set the initialization flag to be true    
	this->Whether_Physical_Group_Info_build=false;           // whether physical group of elements or nodes build

	this->SetNumberOfInputPorts(0);							 // set numer of Input ports to be 0
	this->SetNumberOfOutputPorts(1);						 // set numer of Input ports to be 1

	this->set_VTK_To_ESSI_Elements_Connectivity();           // build VTK_To_ESSI_connectivity map
	this->Build_Meta_Array_Map();                            // build Meta_Dataset_Array_Map

	this->Physical_Element_Group=vtkSmartPointer<vtkDataArraySelection>::New();  // create vtkDataArraySelection object for Physical_Element_Group
	this->Physical_Node_Group=vtkSmartPointer<vtkDataArraySelection>::New();     // create vtkDataArraySelection object for Physical_Node_Group

	// Visualization Ugrid current mesh with array and saved mesh
	Visualization_Current_UGrid_Node_Mesh= vtkSmartPointer<vtkUnstructuredGrid>::New();
	Visualization_Current_UGrid_Gauss_Mesh= vtkSmartPointer<vtkUnstructuredGrid>::New();
	Visualization_UGrid_Node_Mesh= vtkSmartPointer<vtkUnstructuredGrid>::New();
	Visualization_UGrid_Gauss_Mesh= vtkSmartPointer<vtkUnstructuredGrid>::New();

	Build_Inverse_Matrices();								 // build Inverse_Matrices
	Build_Gauss_To_Node_Interpolation_Map();                 // build Gauss_To_Node_Interpolarion_Map
	this->Whether_Writing_Allowed = false;
	this->Whether_Piece_Data_initialized = false;


	id_H5F_CLOSE_STRONG = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fclose_degree(id_H5F_CLOSE_STRONG, H5F_CLOSE_STRONG);

  	id_H5F_READ_ONLY = (H5F_ACC_RDONLY | H5F_ACC_SWMR_READ);
  	id_H5F_READ_WRITE = H5F_ACC_RDWR;


} 

//============================================================================
// RequestData responds to the request made by vtk Pipeleine
//
// * This method is invoked when the time stamp is changed from paraview VCR
//   or something is modified.
// * It returns the current timestamp (in sec) and the reader needs to 
//   return the mesh for that timestamp.
// 
//============================================================================
int pvESSI::RequestData(vtkInformation *vtkNotUsed(request),vtkInformationVector **vtkNotUsed(inputVector),	vtkInformationVector *outputVector){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::RequestData");
#endif

	vtkInformation *Node_Mesh_Info = outputVector->GetInformationObject(0);
	// outInfo->Print(std::cout);

	// Returns the current time of visualization (in sec)
  	Node_Mesh_Current_Time =  Node_Mesh_Info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

{

	vtkInformation *Node_Mesh_Info = outputVector->GetInformationObject(0);

	this->Update_Time_Steps();
	double Time_range[2]={Time[0],Time[Number_of_Time_Steps-1]};

	Node_Mesh_Info->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),Time, this->Number_of_Time_Steps);
	Node_Mesh_Info->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),Time_range,2);

	// this->Superclass::UpdateInformation();
	// Node_Mesh_Info->SetViewTime(Node_Mesh_Current_Time);

	// this->Superclass::Modified();


}


	TimeIndex1=0;
	TimeIndex2=0;

	for (int i=0; i<Number_of_Time_Steps; i++){
		if(Node_Mesh_Current_Time <= Time[i])
		{
			TimeIndex2 = i;
			break;
		}
	}

	if(TimeIndex2==0) TimeIndex1 = 0; else TimeIndex1 = TimeIndex2-1;
	if(TimeIndex2==Number_of_Time_Steps-1) TimeIndex1 = TimeIndex2;

	if(TimeIndex1==TimeIndex2){
		InterpolationFun1 = 0.5;
		InterpolationFun2 = 0.5;
	}
	else{
		InterpolationFun1 = (Time[TimeIndex2]-Node_Mesh_Current_Time)/(Time[TimeIndex2]-Time[TimeIndex1]);
		InterpolationFun2 = (Node_Mesh_Current_Time-Time[TimeIndex1])/(Time[TimeIndex2]-Time[TimeIndex1]);
	}

  	// piece_no := process_id in paralllel visualization
  	// num_pieces := number of process Ids in parallel visualization
	piece_no      = Node_Mesh_Info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
	num_of_pieces = Node_Mesh_Info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

#ifdef DEBUG_MODE
	cout << "piece_no               " << piece_no <<endl;
	cout << "num_of_pieces          " << num_of_pieces << endl;
	cout << "Current_Time[s]        " << Node_Mesh_Current_Time <<endl;
	cout << "TimeIndex1             " << TimeIndex1 << endl;
	cout << "TimeIndex2             " << TimeIndex2 << endl;
	cout << "InterpolationFun1      " << InterpolationFun1 << endl;
	cout << "InterpolationFun2      " << InterpolationFun2 << endl;
#endif

	// ================================================================================
	// Start_Domian_Number = -1, means that the nothing needs to be done
	//       Example:- Sequential or single slave ESSI output Visualization, 
	//
	// Also, single_file_visualization_mode is enabled if process Id 
	//       number of process' used is equal to zero or process id >0
	//
	// In single visualization mode, piece_no  is equal to 0;
	//
	// In paralell or sequential mode when one or more files needs to be 
	//       combined, loop is required to get the desired data. It is readers 
	//		 responsibility to return the full mesh in all cases. 
	//===============================================================================              

	int Start_Domian_Number =-1;
	int End_Domain_Number   = 0;
	if(single_file_visualization_mode){ // If ony one file needs to be considered 
		if(piece_no>0)
			return 1;
	}
	else{ // if multiple files needs to be read and appeneded to get the given piece data
		Start_Domian_Number = piece_no*ceil(((double)Number_of_Processes_Used)/((double)num_of_pieces));
		End_Domain_Number   = Start_Domian_Number +ceil(((double)Number_of_Processes_Used)/((double)num_of_pieces));
		End_Domain_Number   = End_Domain_Number > (Number_of_Processes_Used-1) ? Number_of_Processes_Used-1:End_Domain_Number; 
	}

	// initializes all the data that is required to form that piece 
	// 		Start_Domian_Number -> correponds to the domain or ESSI output process Id to start reading from
	// 		End_Domain_Number   -> corresponds to the domain or ESSI output provess Id to stop reading at
 	Initialize_Piece_data(Start_Domian_Number,End_Domain_Number);


#ifdef DEBUG_MODE
 	cout << "Start_Domain_Number    " << Start_Domian_Number << endl;
 	cout << "End_Domain_Number      " << End_Domain_Number << endl;
#endif

 	cout << "<<<<pvESSI>>>> Piece No " << piece_no << endl;
	for (int i = Start_Domian_Number; i<End_Domain_Number; i++){

		this->domain_no = i;
		Domain_Initializer(domain_no); // also updates the domian_no

		if (!Whether_Node_Mesh_build[domain_no] and (!Show_Gauss_Mesh_Flag or eigen_mode_on or Enable_Displacement_Probing_Flag)){
			this->Get_Node_Mesh(Visualization_UGrid_Node_Mesh);
			Whether_Node_Mesh_build[domain_no] = true;
		}
		else if (!Whether_Gauss_Mesh_build[domain_no] and Show_Gauss_Mesh_Flag){
			this->Get_Gauss_Mesh(Visualization_UGrid_Gauss_Mesh);
			Whether_Gauss_Mesh_build[domain_no] = true;
				// this->Build_Delaunay3D_Gauss_Mesh(Gauss_Mesh);
		}

		if (!Show_Gauss_Mesh_Flag or eigen_mode_on){
			Visualization_Current_UGrid_Node_Mesh->ShallowCopy(Visualization_UGrid_Node_Mesh);
		}
		if (Show_Gauss_Mesh_Flag){
			Visualization_Current_UGrid_Node_Mesh->ShallowCopy(Visualization_UGrid_Gauss_Mesh);
		}
	  	// int Clength = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	  	// double* Csteps = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

#ifdef DEBUG_MODE
	// cout << "Number_of_Nodes           " << Number_of_Nodes           << endl;
	// cout << "Pseudo_Number_of_Nodes    " << Pseudo_Number_of_Nodes    << endl;
	// cout << "Number_of_Gauss Points    " << Number_of_Gauss_Nodes     << endl;
	// cout << "Number_of_Elements        " << Number_of_Elements        << endl;
	// cout << "Pseudo_Number_of_Elements " << Pseudo_Number_of_Elements << endl;
#endif
		if(eigen_mode_on){ // Only if eigen mode visualizatio is on
			Build_Eigen_Modes_Node_Attributes(Visualization_Current_UGrid_Node_Mesh, TimeIndex1, TimeIndex2, InterpolationFun1, InterpolationFun2);	
		}
		else{ // other than that 

			if(!Show_Gauss_Mesh_Flag){ // if gauss mesh is disabled

				// Builds Node Attributes
				Build_Node_Attributes(Visualization_Current_UGrid_Node_Mesh, TimeIndex1, TimeIndex2, InterpolationFun1, InterpolationFun2);

				// Builds Gauss To Node Interpolation if it is enabled
				if(Enable_Gauss_To_Node_Interpolation_Flag) Build_Node_Stress(Visualization_Current_UGrid_Node_Mesh,  TimeIndex1, TimeIndex2, InterpolationFun1, InterpolationFun2);	

				// If Physical Group Selection is enabled
				if(Enable_Physical_Node_Group_Selection_Flag or Enable_Physical_Element_Group_Selection_Flag) Build_Physical_Element_Group_Mesh(Visualization_Current_UGrid_Node_Mesh);
			}

			else if(Show_Gauss_Mesh_Flag){ // if gauss mesh is enabled
				Build_Gauss_Attributes(Visualization_Current_UGrid_Node_Mesh, TimeIndex1, TimeIndex2, InterpolationFun1, InterpolationFun2);
			}
		}
	}

	this->Whether_Node_Mesh_Attributes_Initialized=false;
	this->Whether_Gauss_Mesh_Attributes_Initialized = false;
	this->Whether_Node_Mesh_Stress_Attributes_Initialized = false;

	// vtkSmartPointer<vtkUnstructuredGrid>  Visualization =  vtkSmartPointer<vtkUnstructuredGrid>::New();;

	// if(single_file_visualization_mode){
	// 	Chunk_Node_mesh->ShallowCopy(UGrid_Current_Node_Mesh[0]);
	// }
	// else{
	// 	Merge_Mesh(start,end, Chunk_Node_mesh);
	// }



	// get the ouptut pointer to paraview 
	vtkUnstructuredGrid *Output_Node_Mesh = vtkUnstructuredGrid::SafeDownCast(Node_Mesh_Info->Get(vtkDataObject::DATA_OBJECT()));

#ifdef DEBUG_MODE
	vtkIndent indent;
	Visualization_Current_UGrid_Node_Mesh->PrintSelf(std::cout, indent);
#endif

	Output_Node_Mesh->ShallowCopy(Visualization_Current_UGrid_Node_Mesh); // return the unstrutured mesh

	// static bool updateonce = true;
	// if(updateonce){ 
	// 	cout << " this->Superclass::Modified()" << endl;
	// this->Modified();
	// updateonce = false;
	// }


	return 1;
}

// ****************************************************************************
// * This method is called only once for information about time stamp and extent. 
// * This is the method which is called firt after the default constructor is 
// * initialized
// ****************************************************************************

int pvESSI::RequestInformation( vtkInformation *request, vtkInformationVector **vtkNotUsed(inVec), vtkInformationVector* outVec){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::RequestInformation");
#endif

	this->Update_Time_Steps();


	if(this->Enable_Initialization_Flag){this->Initialize();} this->Enable_Initialization_Flag=false;

	vtkInformation* Node_Mesh = outVec->GetInformationObject(0);

	double Time_range[2]={Time[0],Time[Number_of_Time_Steps-1]};

	Node_Mesh->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),Time, this->Number_of_Time_Steps);
	Node_Mesh->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),Time_range,2);
	Node_Mesh->Set(CAN_HANDLE_PIECE_REQUEST(), 1);

	return 1;
}


void pvESSI::Merge_Mesh(int start, int end, vtkSmartPointer<vtkUnstructuredGrid> Mesh){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Merge_Mesh");
#endif

  vtkSmartPointer<vtkAppendFilter> AppendFilter =  vtkSmartPointer<vtkAppendFilter>::New();

  vtkAppendFilter::GlobalWarningDisplayOff();

  // for (int i = start; i<end; i++){
	 //  AppendFilter->AddInputData(UGrid_Current_Node_Mesh[i]);
	 // // vtkIndent indent;
	 // // delaunay3D->PrintSelf(std::cout, indent); 
  // }

  AppendFilter->Update();
  Mesh->ShallowCopy( AppendFilter->GetOutput());

  return;

}


//============================================================================
// This method creates and process the request from vtk Pipeleine
//============================================================================
int pvESSI::ProcessRequest(vtkInformation  *request_type, vtkInformationVector  **inInfo, vtkInformationVector *outInfo){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::ProcessRequest");
	vtkIndent indent;
 	request_type->PrintSelf(std::cout, indent);
#endif

	return this->Superclass::ProcessRequest(request_type, inInfo, outInfo);
}

//============================================================================
// This method prints about reader plugin i.e itself
//============================================================================
void pvESSI::PrintSelf(ostream& os, vtkIndent indent){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::PrintSelf");
#endif

	this->Superclass::PrintSelf(os,indent);
	os << indent << "File Name: " << (this->FileName ? this->FileName : "(none)") << "\n";
	return;
}

//=============================================================================
// This method creates a map for the connectivity order from ESSI connectivity
// to vtk elements connectivity. The connectivity order is stored in
// ESSI_to_VTK_Element map whose key is the number of connectivity nodes. 
//
// !!!I think a better key is needed so that it can be robust.
//=============================================================================/
void pvESSI::set_VTK_To_ESSI_Elements_Connectivity(){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::set_VTK_To_ESSI_Elements_Connectivity");
#endif

	std::vector<int> connectivity_vector;

	ESSI_to_VTK_Element[1] = VTK_VERTEX;		ESSI_to_VTK_Element[2]  = VTK_LINE;					ESSI_to_VTK_Element[4] = VTK_QUAD;
	ESSI_to_VTK_Element[8] = VTK_HEXAHEDRON;	ESSI_to_VTK_Element[20] = VTK_QUADRATIC_HEXAHEDRON;	ESSI_to_VTK_Element[27] = VTK_TRIQUADRATIC_HEXAHEDRON;

	int Node[1] = {0};						connectivity_vector.assign(Node,Node+1);				ESSI_to_VTK_Connectivity[1] = connectivity_vector;	connectivity_vector.clear();
	int Line[2] = {0,1}; 					connectivity_vector.assign(Line,Line+2); 				ESSI_to_VTK_Connectivity[2] = connectivity_vector;	connectivity_vector.clear();
	int Quadrangle[4] = {0,1,2,3}; 			connectivity_vector.assign(Quadrangle,Quadrangle+4); 	ESSI_to_VTK_Connectivity[4] = connectivity_vector;	connectivity_vector.clear();
	int Hexahedron[8] = {4,5,6,7,0,1,2,3}; 	connectivity_vector.assign(Hexahedron,Hexahedron+8);	ESSI_to_VTK_Connectivity[8] = connectivity_vector;	connectivity_vector.clear();
	int Qudratic_hexahedron[20] = {6,5,4,7,2,1,0,3,13,12,15,14,9,8,11,10,18,17,16,19}; 	connectivity_vector.assign(Qudratic_hexahedron,Qudratic_hexahedron+20);	ESSI_to_VTK_Connectivity[20]= connectivity_vector;	connectivity_vector.clear();
	int TriQudratic_hexahedron[27] = {6,5,4,7,2,1,0,3,13,12,15,14,9,8,11,10,18,17,16,19,23,21,22,24,26,25,20}; 	connectivity_vector.assign(TriQudratic_hexahedron,TriQudratic_hexahedron+27);	ESSI_to_VTK_Connectivity[27]= connectivity_vector;	connectivity_vector.clear();
	return;
}

//===========================================================================
// This Method builds node attributes for the interpolation time step index 
// T1 and T2 and interpolation function InFun1 and InFun2 and pushes it to 
// the <vtkUnstructuredGrid> input vtkobject.
//===========================================================================
void pvESSI::Build_Node_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh,  int Time_Index1, int Time_Index2, float Interpolation_Func1,float Interpolation_Func2) {

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Node_Attributes");
#endif

	////////////////////////////////////////////////// For Reference Displacement  ///////////////////////////////////////////////////////////////////////////
	// Finding out the refernce index for displacement 
	int reference_node_mesh_time = Reference_Displacement_Index_Flag;
	if(Reference_Displacement_Index_Flag>=Number_of_Time_Steps) reference_node_mesh_time=Number_of_Time_Steps;
	else if(Reference_Displacement_Index_Flag<=0) reference_node_mesh_time=0;
	float *Node_Generalized_Displacements,*Reference_Node_Generalized_Displacements;

	///////////////////////////////////////////  Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	
	this->INTERPOLATE_FLOAT_Time_Data_From_2_D_Dataset(id_Generalized_Displacements,Time_Index1, Time_Index2, Interpolation_Func1,Interpolation_Func2,&Node_Generalized_Displacements);
	if(Enable_Relative_Displacement_Flag)
	{
		this->FLOAT_Time_Data_From_2_D_Dataset(id_Generalized_Displacements, reference_node_mesh_time, &Reference_Node_Generalized_Displacements);
	}

	/////////////////////////////////////////////////////////////// DataSets Visulization at Nodes //////////////////////////////////////////////////////////////////////////////////////
 	
 	if(Whether_Node_Mesh_Attributes_Initialized==false)
 	{
	 	Generalized_Displacements = vtkSmartPointer<vtkFloatArray>::New(); 
		this->Set_Meta_Array (Meta_Array_Map["Generalized_Displacements"]);
		Generalized_Displacements->Allocate(Total_Number_of_Nodes*3);
		Generalized_Displacements->SetNumberOfValues(Total_Number_of_Nodes*3);
		Generalized_Displacements->FillComponent(0,0);
		Generalized_Displacements->FillComponent(1,0);
		Generalized_Displacements->FillComponent(2,0);

	 	if(Enable_uPU_Visualization_Flag)
	 	{
		 	Fluid_Displacements = vtkSmartPointer<vtkFloatArray>::New(); 
			this->Set_Meta_Array(Meta_Array_Map["Fluid_Displacements"]);
			Fluid_Displacements->Allocate(Total_Number_of_Nodes*3);
			Fluid_Displacements->SetNumberOfValues(Total_Number_of_Nodes*3);
			Fluid_Displacements->FillComponent(0,0);
			Fluid_Displacements->FillComponent(1,0);
			Fluid_Displacements->FillComponent(2,0);

		 	Pore_Pressure = vtkSmartPointer<vtkFloatArray>::New(); 
			this->Set_Meta_Array (Meta_Array_Map["Pore_Pressure"]);
			Pore_Pressure->Allocate(Total_Number_of_Nodes);
			Pore_Pressure->SetNumberOfValues(Total_Number_of_Nodes);
			Pore_Pressure->FillComponent(0,0);
		}

		if(enable_support_reactions){

			Support_Reactions = vtkSmartPointer<vtkFloatArray>::New(); 
			this->Set_Meta_Array (Meta_Array_Map["Support_Reactions"]);
			Support_Reactions->Allocate(Total_Number_of_Nodes*3);
			Support_Reactions->SetNumberOfValues(Total_Number_of_Nodes*3);


			Support_Reactions->FillComponent(0,0);
			Support_Reactions->FillComponent(1,0);
			Support_Reactions->FillComponent(2,0);
		}

		this->Whether_Node_Mesh_Attributes_Initialized=true;
 	}


 	int index_to_generalized_displacement=0;

	int ESSI_NODE_TAG;
	int numDOF;

	
	if(Enable_Relative_Displacement_Flag==false)
	{
		for (int i = 0; i < this->Domain_Number_of_Nodes[this->domain_no]; i++)
		{

			ESSI_NODE_TAG = this->Domain_Node_Map[this->domain_no][i];
			numDOF        = this->Domain_Number_of_Dofs[this->domain_no][i];

			float disp[3]={
				Node_Generalized_Displacements[index_to_generalized_displacement  ],
				Node_Generalized_Displacements[index_to_generalized_displacement+1],
				Node_Generalized_Displacements[index_to_generalized_displacement+2]
			};


			Generalized_Displacements->InsertTypedTuple(this->ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],disp);

			// ************************************************************** Displacement, Acceleration and Velocity -> Calculation Formulae *********************************************************

			// Displacement(t) = D_t;
			// Acceleration(t) = (1/del.t)^2*{ D_(t+del.t) - 2*D_t + D_(t-del.t)};
			// Acceleration(t) = 1/2/del.t*{ D_(t+del.t) - D_(t-del.t)};

			/***************************************************************** Add Acceleration and Velocity Arrays Here ****************************************/

			// uPU Visualization Mode

		 	if(Enable_uPU_Visualization_Flag){

		 		static float disp[3]; disp[0]=0; disp[1]=0; disp[2]=0;
		 		static float pore_pressure; pore_pressure =0;

		 		if(numDOF==7)
		 		{
			 		disp[0]=Node_Generalized_Displacements[index_to_generalized_displacement+4];
			 		disp[1]=Node_Generalized_Displacements[index_to_generalized_displacement+5];
			 		disp[2]=Node_Generalized_Displacements[index_to_generalized_displacement+6];

			 		pore_pressure = Node_Generalized_Displacements[index_to_generalized_displacement+3];

				}

				Fluid_Displacements->InsertTypedTuple(this->ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],disp);
				Pore_Pressure->InsertValue(this->ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],pore_pressure);
		 	}

		 	index_to_generalized_displacement = index_to_generalized_displacement + numDOF;


		}

	}
	else{

		for (int i = 0; i < this->Domain_Number_of_Nodes[this->domain_no]; i++)
		{
			ESSI_NODE_TAG = this->Domain_Node_Map[this->domain_no][i];
			numDOF        = this->Domain_Number_of_Dofs[this->domain_no][i];

			float disp[3]={
				Node_Generalized_Displacements[index_to_generalized_displacement  ] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement  ],
				Node_Generalized_Displacements[index_to_generalized_displacement+1] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement+1],
				Node_Generalized_Displacements[index_to_generalized_displacement+2] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement+2]
			};

			Generalized_Displacements->InsertTypedTuple(ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],disp);


			// uPU Visualization Mode
			if(Enable_uPU_Visualization_Flag){

				static float disp[3]; disp[0]=0; disp[1]=0; disp[2]=0;
		 		static float pore_pressure; pore_pressure =0;

		 		if(numDOF==7)
		 		{	
		 			disp[0]=Node_Generalized_Displacements[index_to_generalized_displacement+4] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement+4];
		 			disp[1]=Node_Generalized_Displacements[index_to_generalized_displacement+5] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement+5];
		 			disp[2]=Node_Generalized_Displacements[index_to_generalized_displacement+6] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement+6];
					
					pore_pressure = Node_Generalized_Displacements[index_to_generalized_displacement+3] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement+3];
				}
				Fluid_Displacements->InsertTypedTuple(this->ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],disp);
				Pore_Pressure->InsertValue(this->ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],pore_pressure);
		 	}

		 	index_to_generalized_displacement = index_to_generalized_displacement + numDOF;
		}

	}

	Node_Mesh->GetPointData()->AddArray(Generalized_Displacements);


	if(Enable_uPU_Visualization_Flag){	
		Node_Mesh->GetPointData()->AddArray(Fluid_Displacements);
		Node_Mesh->GetPointData()->AddArray(Pore_Pressure);
	}

	/////////////////////////////////////////////////////////////// Support Reactions  //////////////////////////////////////////////////////////////////////////////////////

	if(enable_support_reactions){

		float *Reaction_Forces; INTERPOLATE_FLOAT_Time_Data_From_2_D_Dataset(id_Support_Reactions,Time_Index1, Time_Index2, Interpolation_Func1,Interpolation_Func2,&Reaction_Forces);
		int   *Constrained_Nodes; INT_Time_Data_From_1_D_Dataset(id_Constrained_Nodes, &Constrained_Nodes);
		int   *Constrained_DOFs; INT_Time_Data_From_1_D_Dataset(id_Constrained_DOFs, &Constrained_DOFs);


		for(int i = 0; i<Domain_Number_of_Constrained_Dofs[domain_no]; i++){

			if(Constrained_DOFs[i]<3)
				Support_Reactions->SetComponent(this->ESSI_To_VTK_Node_Map[Constrained_Nodes[i]],Constrained_DOFs[i],Reaction_Forces[i]);

		}

		Node_Mesh->GetPointData()->AddArray(Support_Reactions);
		delete [] Reaction_Forces; Reaction_Forces=NULL;
		delete [] Constrained_Nodes; Constrained_Nodes=NULL;
		delete [] Constrained_DOFs; Constrained_DOFs=NULL;

	}

	delete [] Node_Generalized_Displacements; Node_Generalized_Displacements=NULL;
	if(Enable_Relative_Displacement_Flag==true)	delete [] Reference_Node_Generalized_Displacements; Reference_Node_Generalized_Displacements=NULL;

 	return;
}

//===========================================================================
// This Method builds eigen node attributes for the interpolation time step 
// index T1 and T2 and interpolation function InFun1 and InFun2 and pushes 
// it to the <vtkUnstructuredGrid> input vtkobject.
//===========================================================================
void pvESSI::Build_Eigen_Modes_Node_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh, int Time_Index1, int Time_Index2, float Interpolation_Func1,float Interpolation_Func2){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Eigen_Modes_Node_Attributes");
#endif

	///////////////////////////////////////////  Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	
	float *Node_Generalized_Displacements;this->INTERPOLATE_FLOAT_Time_Data_From_2_D_Dataset(id_Eigen_Modes,Time_Index1, Time_Index2, Interpolation_Func1,Interpolation_Func2,&Node_Generalized_Displacements);
	/////////////////////////////////////////////////////////////// DataSets Visulization at Nodes //////////////////////////////////////////////////////////////////////////////////////
 	
 	if(this->Whether_Node_Mesh_Attributes_Initialized==false)
 	{
	 	Generalized_Displacements = vtkSmartPointer<vtkFloatArray>::New(); 
		this->Set_Meta_Array (Meta_Array_Map["Generalized_Displacements"]);
		Generalized_Displacements->Allocate(Total_Number_of_Nodes);
		Generalized_Displacements->SetNumberOfValues(Total_Number_of_Nodes);
		Generalized_Displacements->FillComponent(0,0);
		Generalized_Displacements->FillComponent(1,0);
		Generalized_Displacements->FillComponent(2,0);


		this->Whether_Node_Mesh_Attributes_Initialized==false;
	}

	int index_to_generalized_displacement=0;
	int ESSI_NODE_TAG;
	

	for (int i = 0; i < this->Domain_Number_of_Nodes[this->domain_no]; i++){

		ESSI_NODE_TAG = this->Domain_Node_Map[this->domain_no][i];

		float disp[3]={
			Node_Generalized_Displacements[index_to_generalized_displacement  ],
			Node_Generalized_Displacements[index_to_generalized_displacement+1],
			Node_Generalized_Displacements[index_to_generalized_displacement+2]
		};

		index_to_generalized_displacement = index_to_generalized_displacement + 3;

		Generalized_Displacements->InsertTypedTuple(this->ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],disp);
	}

	Node_Mesh->GetPointData()->AddArray(Generalized_Displacements);

	delete [] Node_Generalized_Displacements;    Node_Generalized_Displacements=NULL;
 	return;
}


//===========================================================================
// This Method builds gauss node attributes for the interpolation time step 
// index T1 and T2 and interpolation function InFun1 and InFun2 and pushes 
// it to the <vtkUnstructuredGrid> input vtkobject.
//===========================================================================
void pvESSI::Build_Gauss_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Gauss_Mesh, int Time_Index1, int Time_Index2, float Interpolation_Func1,float Interpolation_Func2){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Gauss_Attributes");
#endif

  	if(Enable_Displacement_Probing_Flag){
  		Build_ProbeFilter_Gauss_Mesh(Gauss_Mesh,Time_Index1, Time_Index2, Interpolation_Func1,Interpolation_Func2);
  	}

 	if(Whether_Gauss_Mesh_Attributes_Initialized==false)
 	{
		Elastic_Strain = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Elastic_Strain"]);
		Elastic_Strain->Allocate(Total_Number_of_Gauss_Points*6);
		Elastic_Strain->SetNumberOfValues(Total_Number_of_Gauss_Points*6);
		Elastic_Strain->FillComponent(0,0);Elastic_Strain->FillComponent(3,0);
		Elastic_Strain->FillComponent(1,0);Elastic_Strain->FillComponent(4,0);
		Elastic_Strain->FillComponent(2,0);Elastic_Strain->FillComponent(5,0);


		Plastic_Strain = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Plastic_Strain"]);
		Plastic_Strain->Allocate(Total_Number_of_Gauss_Points*6);
		Plastic_Strain->SetNumberOfValues(Total_Number_of_Gauss_Points*6);
		Plastic_Strain->FillComponent(0,0);Plastic_Strain->FillComponent(3,0);
		Plastic_Strain->FillComponent(1,0);Plastic_Strain->FillComponent(4,0);
		Plastic_Strain->FillComponent(2,0);Plastic_Strain->FillComponent(5,0);

		Stress = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Stress"]);
		Stress->Allocate(Total_Number_of_Gauss_Points*6);
		Stress->SetNumberOfValues(Total_Number_of_Gauss_Points*6);
		Stress->FillComponent(0,0);Stress->FillComponent(3,0);
		Stress->FillComponent(1,0);Stress->FillComponent(4,0);
		Stress->FillComponent(2,0);Stress->FillComponent(5,0);

		q = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["q"]);
		q->Allocate(Total_Number_of_Gauss_Points);
		q->SetNumberOfValues(Total_Number_of_Gauss_Points);
		q->FillComponent(0,0);

		p = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["p"]);
		p->Allocate(Total_Number_of_Gauss_Points);
		p->SetNumberOfValues(Total_Number_of_Gauss_Points);
		p->FillComponent(0,0);

		Plastic_Strain_q = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Plastic_Strain_q"]);
		Plastic_Strain_q->Allocate(Total_Number_of_Gauss_Points);
		Plastic_Strain_q->SetNumberOfValues(Total_Number_of_Gauss_Points);
		Plastic_Strain_q->FillComponent(0,0);

		Plastic_Strain_p = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Plastic_Strain_p"]);
		Plastic_Strain_p->Allocate(Total_Number_of_Gauss_Points);
		Plastic_Strain_p->SetNumberOfValues(Total_Number_of_Gauss_Points);
		Plastic_Strain_p->FillComponent(0,0);


		gauss_no_for_attributes = 0;

		Whether_Gauss_Mesh_Attributes_Initialized=true;
	}

	///////////////////////////////////////////  Gauss Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	
	float *Gauss_Outputs; Gauss_Outputs; INTERPOLATE_FLOAT_Time_Data_From_2_D_Dataset(id_Gauss_Outputs,Time_Index1, Time_Index2, Interpolation_Func1,Interpolation_Func2,&Gauss_Outputs);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// cout << "id_Gauss_Outputs  " << id_Gauss_Outputs << endl;
	//////////////////////////////////// Usefull information about tensor components order ///////////////////////////////////////////////////////////////////////////
	// > Symmetric Tensor the order is 0 1 2 1 3 4 2 4 5 so like in a matrix
	// > 0 1 2
	// > 1 3 4
	// > 2 4 6

	// > Assymetriic Tensor tensor the order is  0 1 2 3 4 5 6 7 8
	// > 0 1 2
	// > 3 4 5
	// > 6 7 8
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 	int ngauss, Gauss_Output_Index, index_to_gauss_output =0;
 	double Var_Plastic_p, Var_Plastic_q, Var_p, Var_q;

	for (int i = 0; i < this->Domain_Number_of_Elements[this->domain_no]; i++){

    	ngauss = NUMBER_OF_GAUSS(Element_Desc_Array[this->Domain_Class_Tags[this->domain_no][i]]);

		for(int j=0; j< ngauss ; j++){

			Gauss_Output_Index = index_to_gauss_output*18;

			float El_Strain_Tuple[6] ={
				Gauss_Outputs[Gauss_Output_Index  ],
				Gauss_Outputs[Gauss_Output_Index+3],
				Gauss_Outputs[Gauss_Output_Index+4],
				Gauss_Outputs[Gauss_Output_Index+1],
				Gauss_Outputs[Gauss_Output_Index+5],
				Gauss_Outputs[Gauss_Output_Index+2]
			};

			float Pl_Strain_Tuple[6] ={
				Gauss_Outputs[Gauss_Output_Index+6],
				Gauss_Outputs[Gauss_Output_Index+9],
				Gauss_Outputs[Gauss_Output_Index+10],
				Gauss_Outputs[Gauss_Output_Index+7],
				Gauss_Outputs[Gauss_Output_Index+11],
				Gauss_Outputs[Gauss_Output_Index+8]
			};

			float Stress_Tuple[6] ={
				Gauss_Outputs[Gauss_Output_Index+12],
				Gauss_Outputs[Gauss_Output_Index+15],
				Gauss_Outputs[Gauss_Output_Index+16],
				Gauss_Outputs[Gauss_Output_Index+13],
				Gauss_Outputs[Gauss_Output_Index+17],
				Gauss_Outputs[Gauss_Output_Index+14]
			};

			// cout << "gauss_no_for_attributes " << gauss_no_for_attributes << "---> " << endl;
			// for (int i = 0; i < 6; ++i)
			// {
			// 	cout << "  " << Stress_Tuple[i] ;
			// }
			// cout << endl;

			Var_Plastic_p   = -1.0*(Gauss_Outputs[Gauss_Output_Index+6]+Gauss_Outputs[Gauss_Output_Index+7]+Gauss_Outputs[Gauss_Output_Index+8]);
			Var_Plastic_q   = sqrt(2.0/9.0* (pow(Gauss_Outputs[Gauss_Output_Index+6]-Gauss_Outputs[Gauss_Output_Index+1],7) +
					   					 pow(Gauss_Outputs[Gauss_Output_Index+7]-Gauss_Outputs[Gauss_Output_Index+2],8) +
					   					 pow(Gauss_Outputs[Gauss_Output_Index+8]-Gauss_Outputs[Gauss_Output_Index+1],6) + 6*
					   					 (	pow(Gauss_Outputs[Gauss_Output_Index+9],2)+
					   					 	pow(Gauss_Outputs[Gauss_Output_Index+10],2)+
					   					 	pow(Gauss_Outputs[Gauss_Output_Index+11],2))
					  					 )
							  );

			Var_p   	   = -1.0/3.0*(Gauss_Outputs[Gauss_Output_Index+12]+Gauss_Outputs[Gauss_Output_Index+13]+Gauss_Outputs[Gauss_Output_Index+14]);
			Var_q 		   = sqrt(1.0/2.0* (pow(Gauss_Outputs[Gauss_Output_Index+12]-Gauss_Outputs[Gauss_Output_Index+13],2) +
									   		pow(Gauss_Outputs[Gauss_Output_Index+13]-Gauss_Outputs[Gauss_Output_Index+14],2) +
									   		pow(Gauss_Outputs[Gauss_Output_Index+14]-Gauss_Outputs[Gauss_Output_Index+13],2) + 6*
									   		(	pow(Gauss_Outputs[Gauss_Output_Index+15],2)+
									   			pow(Gauss_Outputs[Gauss_Output_Index+16],2)+
									   			pow(Gauss_Outputs[Gauss_Output_Index+17],2))
									  		)
							 );

			Elastic_Strain->InsertTypedTuple (gauss_no_for_attributes,El_Strain_Tuple);
			Plastic_Strain->InsertTypedTuple (gauss_no_for_attributes,Pl_Strain_Tuple);
			Stress->InsertTypedTuple (gauss_no_for_attributes,Stress_Tuple);
			Plastic_Strain_p->InsertValue(gauss_no_for_attributes,Var_Plastic_p);
			Plastic_Strain_q->InsertValue(gauss_no_for_attributes,Var_Plastic_q);
			p->InsertValue(gauss_no_for_attributes,Var_p);	
			q->InsertValue(gauss_no_for_attributes,Var_q);		
			index_to_gauss_output+=1;

			gauss_no_for_attributes = gauss_no_for_attributes +1; // global gauss no

		}

	}

	Gauss_Mesh->GetPointData()->AddArray(Stress);
  	Gauss_Mesh->GetPointData()->AddArray(Plastic_Strain);
  	Gauss_Mesh->GetPointData()->AddArray(Elastic_Strain);
	Gauss_Mesh->GetCellData()->AddArray(Plastic_Strain_p);
  	Gauss_Mesh->GetCellData()->AddArray(Plastic_Strain_q);
  	Gauss_Mesh->GetCellData()->AddArray(p);
  	Gauss_Mesh->GetCellData()->AddArray(q);

	delete [] Gauss_Outputs; Gauss_Outputs=NULL;

	return;
}


//===========================================================================
// This Method probes node mesh variables at gauss mesh 
// 		By Default [prob_type =0], probes only generalized displacement
// 	 But the user can choose an option [prob_type = 1] to probe all variables
//===========================================================================
void pvESSI::Build_ProbeFilter_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Probe_Input, int Time_Index1, int Time_Index2, float Interpolation_Func1,float Interpolation_Func2){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_ProbeFilter_Gauss_Mesh");
#endif

	vtkSmartPointer<vtkUnstructuredGrid> Probe_Source = vtkSmartPointer<vtkUnstructuredGrid>::New();
	Probe_Source->ShallowCopy(this->Visualization_UGrid_Node_Mesh);

	////////////////////////////////////// For Debugging ////////////////////////////////
	Probe_Source->GetCellData()->RemoveArray("Material_Tag");
	Probe_Source->GetPointData()->RemoveArray("Node_Tag");
	Probe_Source->GetCellData()->RemoveArray("Element_Tag");
	Probe_Source->GetPointData()->RemoveArray("Boundary_Conditions");
	Probe_Source->GetCellData()->RemoveArray("Class_Tag");
	/////////////////////////////////////////////////////////////////////////////////////

	Build_Node_Attributes(Probe_Source, Time_Index1, Time_Index2, Interpolation_Func1,Interpolation_Func2);

	/************* Initializing Probe filter ******************************************/

	vtkSmartPointer<vtkProbeFilter> ProbeFilter = vtkSmartPointer<vtkProbeFilter>::New();
	ProbeFilter->SetSourceData(Probe_Source);
	ProbeFilter->SetInputData(Probe_Input);
	ProbeFilter->Update();

	Probe_Input->ShallowCopy( ProbeFilter->GetOutput());

	  // vtkSmartPointer<vtkGaussianSplatter> splatter = 
	  // vtkSmartPointer<vtkGaussianSplatter>::New();
	  // splatter->SetInputData(Probe_Input);
	  // splatter->SetSampleDimensions(50,50,50);
	  // splatter->SetRadius(0.5);
	  // splatter->ScalarWarpingOff();
	  // splatter->Update();


	// Probe_Input->ShallowCopy(splatter->GetOutput());

  	return;
}


//===========================================================================
// This method creates a Node Mesh i.e a mesh from the given node data 
// and element data from hdf5 file. The function stores the mesh in the 
// given input <vtkUnstructuredGrid> input object.
//===========================================================================
void pvESSI::Get_Node_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Get_Node_Mesh");
#endif

	//////////////////////////////////////////////////// Reading Mesh Data ///////////////////////////////////////////////////////////////////////////////////////////
	float *Node_Coordinates; FLOAT_Time_Data_From_1_D_Dataset(id_Coordinates,&Node_Coordinates);
	////////////////////////////////////////////////////// Reading Node Attributes /////////////////////////////////////////////////////////////////////////////
	int *Constrained_Nodes; INT_Time_Data_From_1_D_Dataset(id_Constrained_Nodes,&Constrained_Nodes);
	int *Constrained_DOFs; INT_Time_Data_From_1_D_Dataset(id_Constrained_DOFs,&Constrained_DOFs);
	//////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////
	int *Element_Material_Tags; INT_Time_Data_From_1_D_Dataset(id_Material_Tags,&Element_Material_Tags);
	int *Element_Partition;

	/////////////////////////////////////////////////// Building Nodes ///////////////////////////////////////////////////////////////////////////////////////////

	if(this->Whether_Node_Mesh_Array_Initialized==false)
	{
		points = vtkSmartPointer<vtkPoints>::New();
		points->SetNumberOfPoints(Total_Number_of_Nodes);

		Node_Mesh->Allocate(Total_Number_of_Nodes);

		ESSI_To_VTK_Node_Map    = new int[Max_Node_Tag];
		ESSI_To_VTK_Element_Map = new int[Max_Node_Tag];
		VTK_To_ESSI_Node_Map    = new int[Total_Number_of_Nodes];
		VTK_To_ESSI_Element_Map = new int[Total_Number_of_Elements];		


		for(int i=0;i<Max_Node_Tag;i++)
			ESSI_To_VTK_Node_Map[i]=-1;
		for(int i=0;i<Max_Element_Tag;i++)
			ESSI_To_VTK_Element_Map[i]=-1;		

		Node_Tag = vtkSmartPointer<vtkIntArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Node_Tag"]);
		Node_Tag->SetNumberOfTuples(Total_Number_of_Nodes);
		Node_Tag->FillComponent(0,0);

		Material_Tag = vtkSmartPointer<vtkIntArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Material_Tag"]);
		Material_Tag->SetNumberOfTuples(Total_Number_of_Elements);
		Material_Tag->FillComponent(0,0);

		Element_Tag = vtkSmartPointer<vtkIntArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Element_Tag"]);
		Element_Tag->SetNumberOfTuples(Total_Number_of_Elements);
		Element_Tag->FillComponent(0,0);

		Class_Tag = vtkSmartPointer<vtkIntArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Class_Tag"]);
		Class_Tag->SetNumberOfTuples(Total_Number_of_Elements);
		Class_Tag->FillComponent(0,0);

		Boundary_Conditions = vtkSmartPointer<vtkIntArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Boundary_Conditions"]);
		Boundary_Conditions->SetNumberOfTuples(Total_Number_of_Nodes);
		Boundary_Conditions->FillComponent (0, 0);
		Boundary_Conditions->FillComponent (1, 0);
		Boundary_Conditions->FillComponent (2, 0);

		if(Number_of_Processes_Used>1){	

			Partition_Info = vtkSmartPointer<vtkIntArray> ::New();
			this->Set_Meta_Array (Meta_Array_Map["Partition_Info"]);
			Partition_Info->Allocate(Total_Number_of_Elements);	
			Partition_Info->SetNumberOfTuples(Total_Number_of_Elements);
			Partition_Info->FillComponent(0,0);	
		}

		node_no=0;
		element_no=0;
		this->Whether_Node_Mesh_Array_Initialized=true;
	}

	if(Number_of_Processes_Used>1){	

		std::string Source_File = GetSourceFile(this->FileName)+"feioutput";
		hid_t id_Source_File = H5Fopen(Source_File.c_str(), id_H5F_READ_ONLY, id_H5F_CLOSE_STRONG);

		id_Element_Partition = H5Dopen(id_Source_File, "/Model/Elements/Partition", H5P_DEFAULT);
		INT_Time_Data_From_1_D_Dataset(id_Element_Partition,&Element_Partition);
		H5Dclose(id_Element_Partition);

		H5Fclose(id_Source_File);
	}


	int index_to_coordinates=0;
	int ESSI_NODE_TAG; 

	for (int i = 0; i < this->Domain_Number_of_Nodes[this->domain_no]; i++){

		ESSI_NODE_TAG = this->Domain_Node_Map[domain_no][i];

		// index = Index_to_Coordinates[ESSI_NODE_TAG];
		if(ESSI_To_VTK_Node_Map[ESSI_NODE_TAG]==-1)
		{
			Node_Tag->InsertValue(node_no,ESSI_NODE_TAG);
			points->InsertPoint(node_no,Node_Coordinates[index_to_coordinates  ],
									Node_Coordinates[index_to_coordinates+1],
									Node_Coordinates[index_to_coordinates+2]
			);

			ESSI_To_VTK_Node_Map[ESSI_NODE_TAG]=node_no;
			VTK_To_ESSI_Node_Map[node_no]=ESSI_NODE_TAG;

			node_no++;
		}
		index_to_coordinates=index_to_coordinates+3;
	}

	Node_Mesh->SetPoints(points);
	Node_Mesh->GetPointData()->AddArray(Node_Tag);

	/////////////////////////////////////////////// Building Boundary Conditions  ///////////////////////////////////////////////////////////////////////////
	for(int i = 0; i<Domain_Number_of_Constrained_Dofs[domain_no]; i++){

		if(Constrained_DOFs[i]<3){
			Boundary_Conditions->SetComponent(ESSI_To_VTK_Node_Map[Constrained_Nodes[i]],Constrained_DOFs[i],1);
		}

	}
	Node_Mesh->GetPointData()->AddArray(Boundary_Conditions);

	///////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	int connectivity_index =0;
	int nnodes = 0;
	int Cell_Type;
	int ESSI_ELEMENT_TAG; 
	int class_tag = 0;

	for (int i = 0; i < Domain_Number_of_Elements[this->domain_no]; i++){

		ESSI_ELEMENT_TAG = this->Domain_Element_Map[domain_no][i];
		class_tag = Domain_Class_Tags[this->domain_no][i];

		nnodes     = NUMBER_OF_NODES(Element_Desc_Array[class_tag]);  // Number of element nodes

		if(not(Show_Hide_Contact_Flag==true and (class_tag==86  or class_tag==87)) and ESSI_To_VTK_Element_Map[ESSI_ELEMENT_TAG]==-1)
		{

			Material_Tag->InsertValue(element_no,Element_Material_Tags[ESSI_ELEMENT_TAG]);
			Element_Tag->InsertValue(element_no,ESSI_ELEMENT_TAG);
			Class_Tag->InsertValue(element_no,class_tag);

			if(Number_of_Processes_Used>1) {Partition_Info->InsertValue(element_no,Element_Partition[ESSI_ELEMENT_TAG]); }

			vtkIdType Vertices[nnodes];
			Cell_Type = ESSI_to_VTK_Element.find(nnodes)->second;
			std::vector<int> Nodes_Connectivity_Order = ESSI_to_VTK_Connectivity.find(nnodes)->second;

			for(int j=0; j<nnodes ; j++){
				Vertices[j] = ESSI_To_VTK_Node_Map[Domain_Node_Map[this->domain_no][(Domain_Connectivity[this->domain_no][connectivity_index + (Nodes_Connectivity_Order[j])])]];
				// cout << "s " << ESSI_To_VTK_Node_Map[ Element_Connectivity[connectivity_index+Nodes_Connectivity_Order[j]] ] << " " << ESSI_To_VTK_Node_Map[Domain_Node_Map[this->domain_no][(Domain_Connectivity[this->domain_no][connectivity_index + (Nodes_Connectivity_Order[j])])]] << endl;

			}

			Node_Mesh->InsertNextCell(Cell_Type, nnodes, Vertices);

			ESSI_To_VTK_Element_Map[ESSI_ELEMENT_TAG]=element_no;
			VTK_To_ESSI_Element_Map[element_no]=ESSI_ELEMENT_TAG;

			element_no = element_no +1;

		}

		connectivity_index = connectivity_index + nnodes;
	}

	Node_Mesh->GetCellData()->AddArray(Material_Tag);
	Node_Mesh->GetCellData()->AddArray(Element_Tag);
	Node_Mesh->GetCellData()->AddArray(Class_Tag);

	if(Number_of_Processes_Used>1) Node_Mesh->GetCellData()->AddArray(Partition_Info);

	if(Number_of_Processes_Used>1)
	{
		delete [] Element_Partition; Element_Partition=NULL;
	}


	delete [] Node_Coordinates; Node_Coordinates = NULL;
	delete [] Constrained_Nodes; Constrained_Nodes = NULL;
	delete [] Constrained_DOFs; Constrained_DOFs = NULL;
	delete [] Element_Material_Tags; Element_Material_Tags = NULL;

	return;
	
}

//===========================================================================
// Builds a mesh correponding to selected physical group of nodes and elements
// Uses VTKSelection filter to generate the mesh 
// Physical Grou of Nodes nand Elements 
//===========================================================================
void  pvESSI::Build_Physical_Element_Group_Mesh(vtkSmartPointer<vtkUnstructuredGrid> NodeMesh){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Physical_Element_Group_Mesh");
#endif

	int number_of_physical_node_group_arrays = GetNumberOfPhysicalElementGroupArrays();
	int number_of_physical_element_group_arrays = GetNumberOfPhysicalElementGroupArrays();

	vtkSmartPointer<vtkIdTypeArray> Element_Ids = vtkSmartPointer<vtkIdTypeArray>::New();
	Element_Ids->SetNumberOfComponents(1); Element_Ids->SetName("Element_Tag");

	vtkSmartPointer<vtkIdTypeArray> Node_Ids = vtkSmartPointer<vtkIdTypeArray>::New();
	Node_Ids->SetNumberOfComponents(1); Node_Ids->SetName("Node_Tag");	

    hid_t id_Physical_Group_Id;
    int number_of_Element_Ids=0,number_of_Node_Ids=0 ;

	//**********/ Physical Element Groups **********************//
	if(Enable_Physical_Element_Group_Selection_Flag)
	{
		for(int k =0; k<number_of_physical_element_group_arrays;k++)
		{
			// cout << GetPhysicalElementGroupArrayStatus(GetPhysicalElementGroupArrayName(k)) << " " << endl;;
			if( GetPhysicalElementGroupArrayStatus(GetPhysicalElementGroupArrayName(k)) >0){
		
				const char *name = GetPhysicalElementGroupArrayName(k);
				id_Physical_Group_Id =  H5Dopen(id_Physical_Element_Groups,name, H5P_DEFAULT);
				int *Individual_Physical_group; INT_Time_Data_From_1_D_Dataset(id_Physical_Group_Id,&Individual_Physical_group);

				// DataSpace = H5Dget_space(id_Physical_Group_Id);
				// H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
				// int Individual_Physical_group[dims1_out[0]];
				// offset1[0]=0;   	MemSpace = H5Screate_simple(1,dims1_out,NULL); 
				// H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset1,NULL,dims1_out,NULL);
				// H5Dread(id_Physical_Group_Id, H5T_NATIVE_INT, MemSpace, DataSpace, H5P_DEFAULT, Individual_Physical_group); 
				// H5Sclose(MemSpace); status=H5Sclose(DataSpace); 

				H5Dclose(id_Physical_Group_Id);

				for(int l=0; l<dims1_out[0];l++){
					Element_Ids->InsertNextValue(Individual_Physical_group[l]);
					number_of_Element_Ids++;
				}

				delete [] Individual_Physical_group; Individual_Physical_group==NULL;

			}
		}
	}

	if(Enable_Physical_Node_Group_Selection_Flag)
	{
		for(int k =0; k<number_of_physical_node_group_arrays;k++)
		{
			// cout << GetPhysicalNodeGroupArrayStatus(GetPhysicalNodeGroupArrayName(k)) << " " << endl;;
			if(GetPhysicalNodeGroupArrayStatus(GetPhysicalNodeGroupArrayName(k))){

				const char *name = GetPhysicalNodeGroupArrayName(k);
				id_Physical_Group_Id =  H5Dopen(id_Physical_Node_Groups,name, H5P_DEFAULT);
				int *Individual_Physical_group; INT_Time_Data_From_1_D_Dataset(id_Physical_Group_Id,&Individual_Physical_group);

				// DataSpace = H5Dget_space(id_Physical_Group_Id);
				// H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
				// int Individual_Physical_group[dims1_out[0]];
				// offset1[0]=0;   	MemSpace = H5Screate_simple(1,dims1_out,NULL); 
				// H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset1,NULL,dims1_out,NULL);
				// H5Dread(id_Physical_Group_Id, H5T_NATIVE_INT, MemSpace, DataSpace, H5P_DEFAULT, Individual_Physical_group); 
				// H5Sclose(MemSpace); status=H5Sclose(DataSpace); 


				H5Dclose(id_Physical_Group_Id);

				for(int l=0; l<dims1_out[0];l++){
					Node_Ids->InsertNextValue(Individual_Physical_group[l]);
					number_of_Node_Ids++;
				}

				delete [] Individual_Physical_group; Individual_Physical_group==NULL;

			}
		}
	}

	//************************************ Element Selection  **************************************/

    vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
	selectionNode->SetFieldType(vtkSelectionNode::CELL);
	selectionNode->SetContentType(vtkSelectionNode::VALUES);
	selectionNode->SetSelectionList(Element_Ids);
	
	vtkSmartPointer<vtkSelection> selection =   vtkSmartPointer<vtkSelection>::New();
	selection->AddNode(selectionNode);
 
	vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();

	extractSelection->SetInputData(0,NodeMesh);
	extractSelection->SetInputData(1,selection);
	extractSelection->Update();

	//************************************ Node Selection ****************************************/

    vtkSmartPointer<vtkSelectionNode> selectionNode2 = vtkSmartPointer<vtkSelectionNode>::New();
	selectionNode2->SetFieldType(vtkSelectionNode::POINT);
	selectionNode2->SetContentType(vtkSelectionNode::VALUES);
	selectionNode2->SetSelectionList(Node_Ids);
	
	vtkSmartPointer<vtkSelection> selection2 =   vtkSmartPointer<vtkSelection>::New();
	selection2->AddNode(selectionNode2);
 
	vtkSmartPointer<vtkExtractSelection> extractSelection2 = vtkSmartPointer<vtkExtractSelection>::New();

	extractSelection2->SetInputData(0,NodeMesh);
	extractSelection2->SetInputData(1,selection2);
	extractSelection2->Update();


	// vtkIndent indent;
	// extractSelection->PrintSelf(std::cout, indent);

	//*********************************** Merging Mesh *****************************************/

	vtkSmartPointer<vtkAppendFilter> AppendFilter =  vtkSmartPointer<vtkAppendFilter>::New();
	vtkAppendFilter::GlobalWarningDisplayOff();
	if(number_of_Element_Ids>0) AppendFilter->AddInputData(extractSelection->GetOutput());
	if(number_of_Node_Ids>0)    AppendFilter->AddInputData(extractSelection2->GetOutput());
	AppendFilter->Update();

	NodeMesh->ShallowCopy( AppendFilter->GetOutput());

}


//===========================================================================
// This Function builds a gauss mesh from the given gauss mesh coordinates  
// Uses Delaunay3D filter to generate the mesh 
//===========================================================================
void pvESSI::Get_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Gauss_Mesh){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Get_Gauss_Mesh");
#endif

	/////////////////////////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////
	float *Element_Gauss_Point_Coordinates; FLOAT_Time_Data_From_1_D_Dataset(id_Gauss_Point_Coordinates,&Element_Gauss_Point_Coordinates);
	////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	if(this->Whether_Gauss_Mesh_Array_Initialized==false)
	{
		points = vtkSmartPointer<vtkPoints>::New();
		points->SetNumberOfPoints(Total_Number_of_Gauss_Points);

		Gauss_Mesh->Allocate(Total_Number_of_Gauss_Points);

		gauss_no=0;
		this->Whether_Gauss_Mesh_Array_Initialized=true;
	}

	vtkIdType onevertex;

	int gauss_point_no=0;
	int ngauss=0;

	for (int i=0; i < Domain_Number_of_Elements[this->domain_no]; i++){

		ngauss     = NUMBER_OF_GAUSS(Element_Desc_Array[Domain_Class_Tags[this->domain_no][i]]);  // Number of gauss nodes

		for(int j=0; j<ngauss ; j++){
			points->InsertPoint(gauss_no, 
								Element_Gauss_Point_Coordinates[3*gauss_point_no  ],
								Element_Gauss_Point_Coordinates[3*gauss_point_no+1],
								Element_Gauss_Point_Coordinates[3*gauss_point_no+2]
							   );

			/// Adding 1D point to mesh 
			onevertex = gauss_no;
			Gauss_Mesh->InsertNextCell(VTK_VERTEX, 1, &onevertex);

			// updating the gauss point no.
			gauss_point_no++;
			gauss_no = gauss_no +1; // global increment of gauss_no
		}

	}

	Gauss_Mesh->SetPoints(points);

	delete [] Element_Gauss_Point_Coordinates; Element_Gauss_Point_Coordinates = NULL;

	return;
	
}


//===========================================================================
// Generate a tetrahedral mesh from the given input points, and stores the   
// newly generated mesh in the same input vtkUnstructuredGrid object.
//===========================================================================
void pvESSI::Build_Delaunay3D_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> GaussMesh){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Delaunay3D_Gauss_Mesh");
#endif

  vtkSmartPointer<vtkDelaunay3D> delaunay3D = vtkSmartPointer<vtkDelaunay3D>::New();

  vtkDelaunay3D::GlobalWarningDisplayOff();

  delaunay3D->SetInputData (GaussMesh);
  delaunay3D->Update();

  // vtkIndent indent;
  // delaunay3D->PrintSelf(std::cout, indent);

  GaussMesh->ShallowCopy( delaunay3D->GetOutput());

  return;

}


//===========================================================================
// This function generates the metadatay attributes array all at one place   
// I am thinking it to remove it. Was meant to be public accessibe to all but
// currently takes lot of parameters. 
//===========================================================================
void pvESSI::Set_Meta_Array( int Meta_Array_Id ){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Set_Meta_Array");
#endif

	switch (Meta_Array_Id){

		case 0:

			Generalized_Displacements->SetName("Generalized_Displacements");
			Generalized_Displacements->SetNumberOfComponents(3);
			Generalized_Displacements->SetComponentName(0,"UX");
			Generalized_Displacements->SetComponentName(1,"UY");
			Generalized_Displacements->SetComponentName(2,"UZ");
			break;  

		case 1:

			Generalized_Velocity->SetName("Generalized_Velocity");
			Generalized_Velocity->SetNumberOfComponents(3);
			Generalized_Velocity->SetComponentName(0,"VX");
			Generalized_Velocity->SetComponentName(1,"VY");
			Generalized_Velocity->SetComponentName(2,"VZ");
			break;

		case 2:

			Generalized_Acceleration->SetName("Generalized_Acceleration");
			Generalized_Acceleration->SetNumberOfComponents(3);
			Generalized_Acceleration->SetComponentName(0,"AX");
			Generalized_Acceleration->SetComponentName(1,"AY");
			Generalized_Acceleration->SetComponentName(2,"AZ");
			break;

		case 3:

			Elastic_Strain->SetName("Elastic_Strain");
			Elastic_Strain->SetNumberOfComponents(6);
			Elastic_Strain->SetComponentName(0,"eps_xx");
			Elastic_Strain->SetComponentName(1,"eps_xy");
			Elastic_Strain->SetComponentName(2,"eps_xz");
			Elastic_Strain->SetComponentName(3,"eps_yy");
			Elastic_Strain->SetComponentName(4,"eps_yz");
			Elastic_Strain->SetComponentName(5,"eps_zz");
			break;

		case 4:

			Plastic_Strain->SetName("Plastic_Strain");
			Plastic_Strain->SetNumberOfComponents(6);
			Plastic_Strain->SetComponentName(0,"eps_xx");
			Plastic_Strain->SetComponentName(1,"eps_xy");
			Plastic_Strain->SetComponentName(2,"eps_xz");
			Plastic_Strain->SetComponentName(3,"eps_yy");
			Plastic_Strain->SetComponentName(4,"eps_yz");
			Plastic_Strain->SetComponentName(5,"eps_zz");
			break;

		case 5:

			Stress->SetName("Stress");
			Stress->SetNumberOfComponents(6);
			Stress->SetComponentName(0,"sig_xx");
			Stress->SetComponentName(1,"sig_xy");
			Stress->SetComponentName(2,"sig_xz");
			Stress->SetComponentName(3,"sig_yy");
			Stress->SetComponentName(4,"sig_yz");
			Stress->SetComponentName(5,"sig_zz");
			break;

		case 6:
			Material_Tag->SetName("Material_Tag");
			Material_Tag->SetNumberOfComponents(1);
			break;

		case 7:
			Total_Energy->SetName("Total_Energy");
			Total_Energy->SetNumberOfComponents(1);
			break;

		case 8:
			Incremental_Energy->SetName("Incremental_Energy");
			Incremental_Energy->SetNumberOfComponents(1);
			break;

		case 9:
			Node_Tag->SetName("Node_Tag");
			Node_Tag->SetNumberOfComponents(1);
			break;

		case 10:
			Element_Tag->SetName("Element_Tag");
			Element_Tag->SetNumberOfComponents(1);
			break;

		case 11:
			q->SetName("q");
			q->SetNumberOfComponents(1);
			break;

		case 12:
			p->SetName("p");
			p->SetNumberOfComponents(1);
			break;

		case 13:
			Plastic_Strain_q->SetName("$\\epsilon^{pl}_{q}$");
			Plastic_Strain_q->SetNumberOfComponents(1);
			break;

		case 14:
			Plastic_Strain_p->SetName("$\\epsilon^{pl}_{p}$");
			Plastic_Strain_p->SetNumberOfComponents(1);
			break;

		case 15:
			Support_Reactions->SetName("Support_Reactions");
			Support_Reactions->SetNumberOfComponents(3);
			Support_Reactions->SetComponentName(0,"FX");
			Support_Reactions->SetComponentName(1,"FY");
			Support_Reactions->SetComponentName(2,"FZ");
			break;

		case 16:
			Boundary_Conditions->SetName("Boundary_Conditions");
			Boundary_Conditions->SetNumberOfComponents(3);
			Boundary_Conditions->SetComponentName(0,"UX");
			Boundary_Conditions->SetComponentName(1,"UY");
			Boundary_Conditions->SetComponentName(2,"UZ");
			break;

		case 17:
			Class_Tag->SetName("Class_Tag");
			Class_Tag->SetNumberOfComponents(1);
			break;

		case 18:
			Partition_Info->SetName("Partition_Info");
			Partition_Info->SetNumberOfComponents(1);
			break;

		case 19:
			Fluid_Displacements->SetName("Fluid_Displacements");
			Fluid_Displacements->SetNumberOfComponents(3);
			Fluid_Displacements->SetComponentName(0,"UX");
			Fluid_Displacements->SetComponentName(1,"UY");
			Fluid_Displacements->SetComponentName(2,"UZ");
			break;

		case 20:
			Pore_Pressure->SetName("Pore_Pressure");
			Pore_Pressure->SetNumberOfComponents(1);
			break;
	}

	return;

}

//===========================================================================
// Initializes important variables taht would be same for all nodes when 
// run in parallel. It initializes the domain variables of the model.
//===========================================================================
void pvESSI::Initialize(){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Initialize");
#endif

	hid_t datasetId;

	/************************* To make Debug on *******************/
	// this->DebugOn();

	/*************************************************************/

	/*********Initializing some parameters of visualization ******/
	this->Build_Map_Status = 0;
	this->Number_of_Strain_Strain_Info = 22;
	this->single_file_visualization_mode = false;


	/*************** Initialize Data Arrays **********************/
	points = vtkSmartPointer<vtkPoints>::New();

    /***************** File_id **********************************/
    this->id_File = H5Fopen(this->FileName, id_H5F_READ_ONLY, id_H5F_CLOSE_STRONG);;  

    /***************** Model Info *******************************/
	datasetId      = H5Dopen(id_File, "/Number_of_Iterations", H5P_DEFAULT);
	H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Sub_Steps );
	H5Dclose(datasetId);

    datasetId = H5Dopen(id_File, "/Number_of_Processes_Used", H5P_DEFAULT);
	H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Processes_Used);
	H5Dclose(datasetId);

    datasetId           = H5Dopen(id_File, "/Process_Number", H5P_DEFAULT);
	H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Process_Number);
	H5Dclose(datasetId);

	/****************** Element Class Description *************************/
	
	// Finding the Element Class Description for each Element type in Real ESSI Simulator
    datasetId = H5Dopen(id_File, "Model/Elements/Element_Class_Desc", H5P_DEFAULT);
    INT_Time_Data_From_1_D_Dataset(datasetId,&Element_Desc_Array);
    H5Dclose(datasetId);

	// Element Class_DESC_ENCODING
	ELE_TAG_DESC_ENCODING = Element_Desc_Array[0];  

	if(ELE_TAG_DESC_ENCODING==-1)
		ELE_TAG_DESC_ENCODING = 235789180;

	// close the file
	H5Fclose(this->id_File);		

	if(Number_of_Processes_Used==1 || Process_Number>0){
		this->single_file_visualization_mode = true;
		this->Number_of_Processes_Used =1;
	}

	// Whether Domain Data has been build
	Domain_Data_Build_Status       = new bool[Number_of_Processes_Used];
	Domain_Write_Status            = new bool[Number_of_Processes_Used];
	Domain_Read_Status             = new bool[Number_of_Processes_Used];
	Domain_Element_Map_Initialized = new bool[Number_of_Processes_Used];
	Domain_Node_Map_Initialized    = new bool[Number_of_Processes_Used];
	Domain_Basic_Info_Initialized  = new bool[Number_of_Processes_Used];

    /*** Element Data ***/
	Domain_Number_of_Elements           = new int [Number_of_Processes_Used]; 
	Domain_Pseudo_Number_of_Elements    = new int [Number_of_Processes_Used];
	Domain_Number_of_Gauss_Points       = new int [Number_of_Processes_Used];
	Domain_Number_of_Connectivity_Nodes = new int [Number_of_Processes_Used];

    /*** Node Data ***/ 
	Domain_Number_of_Nodes            = new int [Number_of_Processes_Used];
	Domain_Pseudo_Number_of_Nodes     = new int [Number_of_Processes_Used];
	Domain_Number_of_Constrained_Dofs = new int [Number_of_Processes_Used];

  	//// 2-Dimensional Containers [Domain_Number][Array] Parameter 
	Domain_Class_Tags                        = new int*[Number_of_Processes_Used];
	Domain_Connectivity                      = new int*[Number_of_Processes_Used];
	Domain_Element_Map                       = new int*[Number_of_Processes_Used];
	Domain_Node_Map                          = new int*[Number_of_Processes_Used];
	Domain_Inverse_Element_Map               = new int*[Number_of_Processes_Used];
	Domain_Inverse_Node_Map                  = new int*[Number_of_Processes_Used];
	Domain_Number_of_Dofs                    = new int*[Number_of_Processes_Used];
	Domain_Number_of_Elements_Shared         = new int*[Number_of_Processes_Used];
	Domain_Number_of_Gauss_Elements_Shared   = new int*[Number_of_Processes_Used];
	Domain_Number_of_Contact_Elements_Shared = new int*[Number_of_Processes_Used];

  	//// Visualization Parameters
	Whether_Node_Mesh_build                  = new bool[Number_of_Processes_Used];
	Whether_Gauss_Mesh_build                 = new bool[Number_of_Processes_Used];
	Domain_Building_of_Map_Status            = new bool[Number_of_Processes_Used];

	for(int i=0; i<Number_of_Processes_Used;i++)
	{
		Whether_Node_Mesh_build[i]=false;
		Whether_Gauss_Mesh_build[i]=false;
		Domain_Data_Build_Status[i]=false;
		Domain_Write_Status[i]=false;
		Domain_Read_Status[i]=false;
		Domain_Node_Map_Initialized[i]=false;
		Domain_Element_Map_Initialized[i]=false;
		Domain_Basic_Info_Initialized[i]=false;
	}

}

//===========================================================================
// Set the building of map status to be true on 
// request from GUI
//===========================================================================
void pvESSI::Set_Domain_Building_of_Map_Status()
{
	for(int i=0; i<Number_of_Processes_Used; i++)
		Domain_Building_of_Map_Status[i]=true;
	this->Modified();
}

//===========================================================================
// Initializes important variables for the given piece visualization
// One piece can contain one or mode domains or process ids.
//===========================================================================
void pvESSI::Initialize_Piece_data(int start, int end)
{

    if(this->Whether_Piece_Data_initialized==false)
    {
#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Initialize_Piece_data");
#endif

	this->Total_Number_of_Elements =0;
	this->Total_Number_of_Nodes =0;
	this->Total_Number_of_Gauss_Points =0;

	this->Max_Node_Tag=0;
	this->Max_Element_Tag=0;

	int NumElements, NumNodes, NumGauss, MaxNode, MaxElement;
	std::string filename;
	hid_t id_Domain_File;

	hid_t datasetId;

	for (int File_No = start; File_No<end; File_No++){

		if(File_No>=0){
			std::string Source_File = GetSourceFile(this->FileName);
			std::stringstream ss;
			int digits = Number_of_Processes_Used > 0 ? (int) log10 ((double) Number_of_Processes_Used) + 1 : 1;
			ss << setfill('0') << setw(digits) << File_No+1;
			filename = Source_File + ss.str()+".feioutput";
		}
		else{
			filename = this->FileName;
		}

		id_Domain_File = H5Fopen(filename.c_str(), id_H5F_READ_ONLY, id_H5F_CLOSE_STRONG);

		datasetId   = H5Dopen(id_Domain_File, "/Number_of_Elements", H5P_DEFAULT);
		H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&NumElements );
		H5Dclose(datasetId);
		this->Total_Number_of_Elements = this->Total_Number_of_Elements + NumElements;

		datasetId   = H5Dopen(id_Domain_File, "/Number_of_Nodes", H5P_DEFAULT);
		H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&NumNodes );
		H5Dclose(datasetId);
		this->Total_Number_of_Nodes = this->Total_Number_of_Nodes + NumNodes;

		datasetId   = H5Dopen(id_Domain_File, "/Number_of_Gauss_Points", H5P_DEFAULT);
		H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&NumGauss );
		H5Dclose(datasetId);
		this->Total_Number_of_Gauss_Points = this->Total_Number_of_Gauss_Points + NumGauss;
		
		datasetId = H5Dopen(id_Domain_File, "/Model/Elements/Index_to_Connectivity", H5P_DEFAULT);
		DataSpace = H5Dget_space(datasetId);
		H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
		MaxElement = dims1_out[0];
		H5Sclose(DataSpace); H5Dclose(datasetId);
		if(this->Max_Element_Tag<MaxElement)this->Max_Element_Tag=MaxElement;

		datasetId = H5Dopen(id_Domain_File, "Model/Nodes/Number_of_DOFs", H5P_DEFAULT);
		DataSpace = H5Dget_space(datasetId);
		H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
		MaxNode = dims1_out[0];
		H5Sclose(DataSpace); H5Dclose(datasetId);
		if(this->Max_Node_Tag<MaxNode)this->Max_Node_Tag=MaxNode;

		H5Fclose(id_Domain_File);
	}

#ifdef DEBUG_MODE
	cout << "Total_Number_of_Elements           " << Total_Number_of_Elements << endl;
	cout << "Total_Number_of_Nodes              " << Total_Number_of_Nodes  << endl;
	cout << "Total_Number_of_Gauss_Points       " << Total_Number_of_Gauss_Points << endl;
	cout << "Max_Node_Tag                       " << Max_Node_Tag  << endl;
	cout << "Max_Element_Tag                    " << Max_Element_Tag << endl;
#endif

	this->Whether_Piece_Data_initialized=true;
	}

}


//===========================================================================
// Initializes important variables for each domain or process Id of 
// Real-ESSI outout
//===========================================================================

void pvESSI::Domain_Initializer(int  Domain_Number){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Domain_Initializer");
#endif

	if(Domain_Number>=0){
		std::string Source_File = GetSourceFile(this->FileName);
		std::stringstream ss;
		int digits = Number_of_Processes_Used > 0 ? (int) log10 ((double) Number_of_Processes_Used) + 1 : 1;
		ss << setfill('0') << setw(digits) << Domain_Number+1;
		Domain_FileName = Source_File + ss.str()+".feioutput";
	}
	else{
		Domain_FileName = this->FileName; 
		domain_no =0;               // setting the index to be zero for Single_Domain_Visualization_Mode
	}

	const char* filename = Domain_FileName.c_str();

	if(!Whether_Physical_Group_Info_build)
	{
		//********************* Building Avialable Physical Groups ***************************************************************/
		// getting the master file containing physical group
		std::string Source_File = GetSourceFile(this->FileName)+"feioutput";
		hid_t id_Source_File = H5Fopen(Source_File.c_str(), id_H5F_READ_ONLY, id_H5F_CLOSE_STRONG);

	    this->id_Physical_Element_Groups = H5Gopen(id_Source_File, "/Model/Physical_Groups/Physical_Element_Groups", H5P_DEFAULT);
	    this->id_Physical_Node_Groups    = H5Gopen(id_Source_File, "/Model/Physical_Groups/Physical_Node_Groups", H5P_DEFAULT);

	    status = H5Ovisit (id_Physical_Element_Groups, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, NULL);
		for(int i=0; i<Physical_Group_Container.size(); ++i){ const char *name = Physical_Group_Container[i].c_str(); this->Physical_Element_Group->AddArray(name);}
		if(!Enable_Physical_Element_Group_Selection_Flag){this->Physical_Element_Group->DisableAllArrays();} Physical_Group_Container.clear();

	    status = H5Ovisit (id_Physical_Node_Groups, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, NULL);
	    for(int i=0; i<Physical_Group_Container.size(); ++i){ const char *name = Physical_Group_Container[i].c_str(); this->Physical_Node_Group->AddArray(name);}
	    if(!Enable_Physical_Node_Group_Selection_Flag){this->Physical_Node_Group->DisableAllArrays();} Physical_Group_Container.clear();

	    Whether_Physical_Group_Info_build = true;;

		H5Gclose(id_Physical_Element_Groups); 
		H5Gclose(id_Physical_Node_Groups);
		H5Fclose(id_Source_File); 
	}

	/****************************************************************************************
	* We want to open the file and keep it open until we have read all the data what we want
	* So lets go wohoo!!!
	****************************************************************************************/

	/***************** File_id **********************************/
	cout << "<<<<pvESSI>>>>  " << filename << " Time_Step " << Node_Mesh_Current_Time << endl << endl;;
	// cout << "Domain_Building_of_Map_Status[domain_no] " << Domain_Building_of_Map_Status[domain_no] << endl;
	
	H5Fclose(this->id_File); // close the previous file
	this->id_File = H5Fopen(filename, id_H5F_READ_ONLY, id_H5F_CLOSE_STRONG);

	Domain_Write_Status[domain_no] = true;

	if(this->Domain_Data_Build_Status[domain_no]==false or Domain_Building_of_Map_Status[domain_no]==true)
	{

		H5Eset_auto (NULL, NULL, NULL);  // To stop HDF% from printing error message
		this->id_pvESSI = H5Gopen(id_File, "/pvESSI", H5P_DEFAULT);  

		if(Domain_Building_of_Map_Status[domain_no]==true)
		{
			this->Build_Local_Domain_Maps(domain_no);
			this->Write_Local_Domain_Maps(domain_no);
			Domain_Building_of_Map_Status[domain_no]=false;
		}
		else{
			if(id_pvESSI>0) {// it can be read, so lets go and read it
				cout << "<<<<pvESSI>>>>  Maps are build HURRAY!!! \n" << endl;
				this->Read_Local_Domain_Maps(domain_no);
				this->Domain_Read_Status[domain_no] = true;
			}
			else{// We need to build the pvESSI folder
				this->Build_Local_Domain_Maps(domain_no);
				this->Write_Local_Domain_Maps(domain_no);

			}
		}

		H5Gclose(this->id_pvESSI);
		this->Domain_Data_Build_Status[domain_no]=true;
	}

#ifdef DEBUG_MODE
	cout << "Domain_Number::                 " << domain_no << endl;
	cout << "========================== Before =========================" << endl;
	cout << "Domain_Data_Build_Status        " << Domain_Data_Build_Status[domain_no] << endl;
	cout << "Domain_Write_Status             " << Domain_Write_Status[domain_no] << endl;
	cout << "Domain_Read_Status              " << Domain_Read_Status[domain_no]  << endl;
	cout << "Domain_Building_of_Map_Status[domain_no]    " << Domain_Building_of_Map_Status[domain_no] << endl;
#endif

	this->openDatasetIds();
}


void pvESSI::openDatasetIds(){

	if(eigen_mode_on){
		this->id_Eigen_Frequencies					   = H5Dopen(id_File,"Eigen_Mode_Analysis/Eigen_Frequencies",H5P_DEFAULT); 	
		this->id_Eigen_Modes 						   = H5Dopen(id_File,"Eigen_Mode_Analysis/Eigen_Modes",H5P_DEFAULT);
	}
	/************************** Open Dataset for reading and *******************************************************************/
	this->id_Constrained_Nodes                 = H5Dopen(id_File, "Model/Nodes/Constrained_Nodes", H5P_DEFAULT);
	this->id_Material_Tags                     = H5Dopen(id_File, "Model/Elements/Material_Tags", H5P_DEFAULT);
	this->id_Physical_Element_Groups           = H5Gopen(id_File, "Model/Physical_Groups/Physical_Element_Groups", H5P_DEFAULT);
	this->id_Physical_Node_Groups              = H5Gopen(id_File, "Model/Physical_Groups/Physical_Node_Groups", H5P_DEFAULT);
	  
	// /**************** Element Info ******************************/
	this->id_Gauss_Point_Coordinates = H5Dopen(id_File, "Model/Elements/Gauss_Point_Coordinates", H5P_DEFAULT);
	// this->id_Element_Outputs = H5Dopen(id_File, "Model/Elements/Element_Outputs", H5P_DEFAULT);
	this->id_Gauss_Outputs = H5Dopen(id_File, "Model/Elements/Gauss_Outputs", H5P_DEFAULT);
	// this->id_Substep_Outputs = H5Dopen(id_File, "Model/Elements/Substep_Outputs", H5P_DEFAULT);

	/**************** Node Info ******************************/
	this->id_Constrained_DOFs = H5Dopen(id_File, "Model/Nodes/Constrained_DOFs", H5P_DEFAULT);
	this->id_Coordinates = H5Dopen(id_File, "Model/Nodes/Coordinates", H5P_DEFAULT);
	this->id_Index_to_Coordinates = H5Dopen(id_File, "/Model/Nodes/Index_to_Coordinates", H5P_DEFAULT);
	this->id_Generalized_Displacements = H5Dopen(id_File, "Model/Nodes/Generalized_Displacements", H5P_DEFAULT);
	this->id_Index_to_Displacements = H5Dopen(id_File, "Model/Nodes/Index_to_Generalized_Displacements", H5P_DEFAULT);
	this->id_Support_Reactions = H5Dopen(id_File,"/Model/Nodes/Support_Reactions", H5P_DEFAULT);
	if(id_Support_Reactions>0) enable_support_reactions=true; else enable_support_reactions=false;

}

void pvESSI::closeDatasetIds(){

	if(eigen_mode_on){
		H5Dclose(this->id_Eigen_Frequencies);
		H5Dclose(this->id_Eigen_Modes);
	}

	/************************** Open Dataset for reading and *******************************************************************/
	  H5Dclose(this->id_Constrained_Nodes);
	  H5Dclose(this->id_Material_Tags);
	  H5Dclose(this->id_Physical_Element_Groups);
	  H5Dclose(this->id_Physical_Node_Groups);
	  
	  // /**************** Element Info ******************************/
	  H5Dclose(this->id_Gauss_Point_Coordinates);
	  // H5Dclose(this->id_Element_Outputs);
	  H5Dclose(this->id_Gauss_Outputs);
	  // this->id_Substep_Outputs = H5Dopen(id_File, "Model/Elements/Substep_Outputs", H5P_DEFAULT);

	  /**************** Node Info ******************************/
	  H5Dclose(this->id_Constrained_DOFs);
	  H5Dclose(this->id_Coordinates);
	  H5Dclose(this->id_Index_to_Coordinates);
	  H5Dclose(this->id_Generalized_Displacements);
	  H5Dclose(this->id_Index_to_Displacements);
	  H5Dclose(this->id_Support_Reactions);
}



//??===========================================================================
//?? Initializes important variables taht would be same for all nodes when 
//?? run in parallel. It initializes the domain variables of the model.
//??===========================================================================

void pvESSI::Close_File(){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Close_File");
#endif

	H5Fclose(this->id_File);  
}


/*****************************************************************************
* Builds maps in hdf5 output file for efficient visualization
* Creates a maps group in the files with the following data set 
*	Node_Map: From global/(input) node number to reduced node numbers
*	Element_Map : From global(input) element number to reduced element numbers
* 	Inverse_Node_Map: From reduced node number to global(input) node number 
* 	Inverse_Element_Map: From reduced element number to global(input) element number 
*****************************************************************************/

    // status = H5Ovisit (id_Physical_Element_Groups, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, NULL);
    // for(int i=0; i<Physical_Group_Container.size(); ++i){const char *name = Physical_Group_Container[i].c_str(); ; this->Physical_Element_Group->AddArray(name);}

    // status = H5Ovisit (id_Physical_Node_Groups, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, NULL);
    // for(int i=0; i<Physical_Group_Container.size(); ++i){const char *name = Physical_Group_Container[i].c_str(); this->Physical_Node_Group->AddArray(name);}


//##########################################################################################
/* Build Local pvESSI folder data for each domain and store it in the memory 
/* the write function would write it back once the file is writable.
/* Needed for real time visualization 
/* This function should be called only once for each domain when pvESSI folder is not written
/* If pvESSI folder is written, there is no need to call this function or only 
/* on special request must be called.                                               */
/* it should be independendent does not depend upon any variable */
void pvESSI::Build_Local_Domain_Maps(int domain_no){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Local_Domain_Maps");
#endif

	cout << "<<<<pvESSI>>>>  Maps Not Build, So lets go and build it first \n " << endl;


	hid_t datasetId;

	// Finiding Pseudo number of nodes 
   	datasetId = H5Dopen(id_File, "Model/Nodes/Number_of_DOFs", H5P_DEFAULT);
    DataSpace = H5Dget_space(datasetId);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
	this->Domain_Pseudo_Number_of_Nodes[domain_no] = dims1_out[0];
	H5Sclose(DataSpace);
	H5Dclose(datasetId);

	// Finiding Pseudo number of elements 
	datasetId = H5Dopen(id_File, "Model/Elements/Class_Tags", H5P_DEFAULT);
	DataSpace = H5Dget_space(datasetId);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
	this->Domain_Pseudo_Number_of_Elements[domain_no] = dims1_out[0];
	H5Sclose(DataSpace);
	H5Dclose(datasetId);

	// Finiding Number of connectivity nodes
	datasetId = H5Dopen(id_File, "Model/Elements/Connectivity", H5P_DEFAULT);
	DataSpace = H5Dget_space(datasetId);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
	this->Domain_Number_of_Connectivity_Nodes[domain_no] = dims1_out[0];
	H5Sclose(DataSpace);
	H5Dclose(datasetId);

	// Finding Number of Constrained nodes
	datasetId = H5Dopen(id_File, "/Model/Nodes/Constrained_DOFs", H5P_DEFAULT);
	DataSpace = H5Dget_space(datasetId);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
	this->Domain_Number_of_Constrained_Dofs[domain_no] = dims1_out[0];
	H5Sclose(DataSpace);
	H5Dclose(datasetId);

	// Finiding total Number of gauss Points
	datasetId   = H5Dopen(id_File, "/Number_of_Gauss_Points", H5P_DEFAULT);
	H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Gauss_Points );
	this->Domain_Number_of_Gauss_Points[domain_no]=this->Number_of_Gauss_Points;
	H5Dclose(datasetId);

	// Finiding total Number of elements
	datasetId   = H5Dopen(id_File, "/Number_of_Elements", H5P_DEFAULT);
	H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Elements );
	this->Domain_Number_of_Elements[domain_no]=this->Number_of_Elements;
	H5Dclose(datasetId);

	// Finiding total Number of nodes
	datasetId  = H5Dopen(id_File, "/Number_of_Nodes", H5P_DEFAULT);
	H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Nodes );
	this->Domain_Number_of_Nodes[domain_no]=this->Number_of_Nodes;
	H5Dclose(datasetId);

#ifdef DEBUG_MODE
	cout << "this->Domain_Pseudo_Number_of_Nodes["<<domain_no<<"]       " << this->Domain_Pseudo_Number_of_Nodes[domain_no] << endl; 
	cout << "this->Domain_Pseudo_Number_of_Elements["<<domain_no<<"]    " << this->Domain_Pseudo_Number_of_Elements[domain_no] << endl; 
	cout << "this->Domain_Number_of_Connectivity_Nodes["<<domain_no<<"] " << this->Domain_Number_of_Connectivity_Nodes[domain_no] << endl; 
	cout << "this->Domain_Number_of_Gauss_Points["<<domain_no<<"]       " << this->Domain_Number_of_Gauss_Points[domain_no] << endl; 
	cout << "this->Domain_Number_of_Elements["<<domain_no<<"]           " << this->Domain_Number_of_Elements[domain_no] << endl; 
	cout << "this->Domain_Number_of_Nodes["<<domain_no<<"]              " << this->Domain_Number_of_Nodes[domain_no] << endl; 
#endif

	//************** Building Node Map datset *******************************//
	this->Domain_Node_Map[domain_no]         = new int[this->Domain_Number_of_Nodes[domain_no]];
	this->Domain_Inverse_Node_Map[domain_no] = new int[this->Domain_Pseudo_Number_of_Nodes[domain_no]];
	this->Domain_Number_of_Dofs[domain_no]   = new int[this->Domain_Number_of_Nodes[domain_no]];

	datasetId = H5Dopen(id_File, "Model/Nodes/Number_of_DOFs", H5P_DEFAULT);
   	int *Number_of_DOFs; INT_Time_Data_From_1_D_Dataset(datasetId, &Number_of_DOFs);
   	H5Dclose(datasetId);

	int index   = -1;
	int Node_Id = 0;
	while(++index < this->Domain_Pseudo_Number_of_Nodes[domain_no]){

	     if(Number_of_DOFs[index]!=-1){

	     	this->Domain_Node_Map[domain_no][Node_Id] = index;
	     	this->Domain_Number_of_Dofs[domain_no][Node_Id] = Number_of_DOFs[index];
	     	this->Domain_Inverse_Node_Map[domain_no][index] = Node_Id;
	     	Node_Id = Node_Id +1; // increment node Id by 1 
	     }
	     else{

	     	this->Domain_Inverse_Node_Map[domain_no][index]=-1; // means element does not exist
	     }
	}

	//************** Building Element Map datset *******************************//
	this->Domain_Element_Map[domain_no]                     = new int[Domain_Number_of_Elements[domain_no]];
	this->Domain_Connectivity[domain_no]                    = new int[Domain_Number_of_Connectivity_Nodes[domain_no]];
	this->Domain_Class_Tags[domain_no]                      = new int[Domain_Number_of_Elements[domain_no]];
	this->Domain_Inverse_Element_Map[domain_no]             = new int[Domain_Pseudo_Number_of_Elements[domain_no]];
	this->Domain_Number_of_Elements_Shared[domain_no]       = new int[Domain_Number_of_Nodes[domain_no]];
	this->Domain_Number_of_Gauss_Elements_Shared[domain_no] = new int[Domain_Number_of_Nodes[domain_no]];

	for (int i = 0 ; i<Domain_Number_of_Nodes[domain_no]; i++){
		this->Domain_Number_of_Elements_Shared[domain_no][i]       = 0;
		this->Domain_Number_of_Gauss_Elements_Shared[domain_no][i] = 0;
	}

	datasetId = H5Dopen(id_File, "Model/Elements/Class_Tags", H5P_DEFAULT);
	int *Element_Class_Tags;   INT_Time_Data_From_1_D_Dataset(datasetId, &Element_Class_Tags);
	H5Dclose(datasetId);

	datasetId = H5Dopen(id_File, "Model/Elements/Connectivity", H5P_DEFAULT);
	int *Element_Connectivity; INT_Time_Data_From_1_D_Dataset(datasetId, &Element_Connectivity);
	H5Dclose(datasetId);

	std::map<int,double**>::iterator it;

	int nnodes     = 0;
	int class_tag  = 0;
	int class_desc = 0;
	int ngauss     = 0;
	int index_to_connectivity = 0 ;

	int index_2;
	index   = -1;
	int Element_Id = 0;

	while(++index < this->Domain_Pseudo_Number_of_Elements[domain_no]){

		class_tag  = Element_Class_Tags[index];
		class_desc = Element_Desc_Array[class_tag];
        nnodes     = NUMBER_OF_NODES(class_desc); // Number of element nodes
        ngauss     = NUMBER_OF_GAUSS(class_desc); // Number of gauss nodes

	     if(class_tag!=-1){

	     	this->Domain_Element_Map[domain_no][Element_Id] = index;
	     	this->Domain_Class_Tags[domain_no][Element_Id]  = class_tag;
	     	this->Domain_Inverse_Element_Map[domain_no][index] = Element_Id;
	     	Element_Id = Element_Id +1; 


	     	// cout << " this->Domain_Element_Map[domain_no][Element_Id] " << index <<endl; 
	     	// cout << " this->Domain_Inverse_Element_Map[domain_no][index] " << Element_Id << endl;
	     	// cout << " this->Domain_Class_Tags[domain_no][Element_Id] " << class_tag << endl << endl;;

	     	if(ngauss>0){

	     		it = Gauss_To_Node_Interpolation_Map.find(class_tag);

				   if(it != Gauss_To_Node_Interpolation_Map.end()){

				   		index_2= -1;
						while(++index_2 <nnodes){

							Node_Id = this->Domain_Inverse_Node_Map[domain_no][Element_Connectivity[index_2+index_to_connectivity]];
							this->Domain_Number_of_Gauss_Elements_Shared[domain_no][Node_Id] = this->Domain_Number_of_Gauss_Elements_Shared[domain_no][Node_Id]+1;
							this->Domain_Number_of_Elements_Shared[domain_no][Node_Id]      = this->Domain_Number_of_Elements_Shared[domain_no][Node_Id]+1;
							this->Domain_Connectivity[domain_no][index_2+index_to_connectivity] = Node_Id;
						}
				   }
				   else{
				   	   	
				   	   	index_2= -1;
						while(++index_2 <nnodes){
							Node_Id = this->Domain_Inverse_Node_Map[domain_no][Element_Connectivity[index_2+index_to_connectivity]];
							this->Domain_Number_of_Elements_Shared[domain_no][Node_Id]      = this->Domain_Number_of_Elements_Shared[domain_no][Node_Id]+1;
							this->Domain_Connectivity[domain_no][index_2+index_to_connectivity] = Node_Id;
						}

				   	   cout << "<<<<[pvESSI]>>>> Build_Maps:: Warning!! Gauss to Node Interpolation not implemented for element of Class_Tag  " << Element_Class_Tags[index] << endl;

				   }

	     	}
		     else{

			   	   	index_2= -1;
					while(++index_2 <nnodes){
						Node_Id = this->Domain_Inverse_Node_Map[domain_no][Element_Connectivity[index_2+index_to_connectivity]];
						this->Domain_Number_of_Elements_Shared[domain_no][Node_Id]      = this->Domain_Number_of_Elements_Shared[domain_no][Node_Id]+1;
						this->Domain_Connectivity[domain_no][index_2+index_to_connectivity] = Node_Id;
					}

		     }

		    index_to_connectivity = index_to_connectivity + nnodes;
	     }
	     else{

	     	this->Domain_Inverse_Element_Map[domain_no][index]=-1;
	     }
	}

	delete[] Element_Class_Tags; Element_Class_Tags = NULL;
	delete[] Element_Connectivity; Element_Connectivity = NULL;
	delete[] Number_of_DOFs; Number_of_DOFs = NULL;

	// Set Domain Data Build Status to True 
	this->Domain_Data_Build_Status[domain_no]=true;

	cout << "<<<<pvESSI>>>>  Maps are now built \n " << endl;

}


//######################################################################
// Write the Local Domain Maps to file if it has already been build 
// inside pvESSI folder
void pvESSI::Write_Local_Domain_Maps(int domain_no){


#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Write_Local_Domain_Maps");
#endif

	if(this->Domain_Write_Status[domain_no]==false){
		return;
	}

	// Check if Local_Domains Maps were build else build it again //
	if(!Domain_Data_Build_Status[domain_no]) 
		Build_Local_Domain_Maps(domain_no);

	// Close the file and reopen in write mode and again close te file 
	// and reopen in read only mode
	H5Fclose(this->id_File);
	const char * filename = Domain_FileName.c_str();
	this->id_File = H5Fopen(filename, id_H5F_READ_WRITE, id_H5F_CLOSE_STRONG);

	if(this->id_File )
	{	
		cout << "<<<pvESSI>>>  Writing Maps to File " << endl;


		int ret, status; 
	    ret = H5Ldelete(this->id_File, "/pvESSI", H5P_DEFAULT);
	    status = H5Lexists(this->id_File, "/pvESSI", H5P_DEFAULT);

		hid_t datasetId;

		/*********************************************** First Need to create Folders *********************************************************************/
		id_pvESSI = H5Gcreate(id_File, "/pvESSI", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
		H5Gclose(id_pvESSI); 

		datasetId = H5Gcreate(id_File, "/pvESSI/Field_at_Nodes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
		H5Gclose(datasetId); 

		dims1_out[0]= this->Domain_Number_of_Nodes[domain_no];
		DataSpace = H5Screate_simple(1, dims1_out, NULL);

		datasetId= H5Dcreate(id_File,"pvESSI/Node_Map",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(DataSpace);
		H5Dwrite(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,this->Domain_Node_Map[domain_no]);
		H5Dclose(datasetId);

		DataSpace = H5Screate_simple(1, dims1_out, NULL);
		datasetId = H5Dcreate(id_File,"pvESSI/Number_of_DOfs",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(DataSpace);
		H5Dwrite(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Domain_Number_of_Dofs[domain_no]);
		H5Dclose(datasetId);

		DataSpace = H5Screate_simple(1, dims1_out, NULL);
		datasetId = H5Dcreate(id_File,"pvESSI/Number_of_Elements_Shared",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(DataSpace);
		H5Dwrite(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Domain_Number_of_Elements_Shared[domain_no]);
		H5Dclose(datasetId);

		DataSpace = H5Screate_simple(1, dims1_out, NULL);
		datasetId = H5Dcreate(id_File,"pvESSI/Number_of_Gauss_Elements_Shared",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(DataSpace);
		H5Dwrite(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Domain_Number_of_Gauss_Elements_Shared[domain_no]);
		H5Dclose(datasetId);

		dims1_out[0]= this->Domain_Number_of_Elements[domain_no]; 
		DataSpace = H5Screate_simple(1, dims1_out, NULL);
		datasetId = H5Dcreate(id_File,"pvESSI/Element_Map",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(DataSpace);
		H5Dwrite(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Domain_Element_Map[domain_no]);
		H5Dclose(datasetId);

		DataSpace = H5Screate_simple(1, dims1_out, NULL);
		datasetId = H5Dcreate(id_File,"pvESSI/Class_Tags",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(DataSpace);
		H5Dwrite(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Domain_Class_Tags[domain_no]);
		H5Dclose(datasetId);

		dims1_out[0]= this->Domain_Number_of_Connectivity_Nodes[domain_no]; 
		DataSpace = H5Screate_simple(1, dims1_out, NULL);
		datasetId = H5Dcreate(id_File,"pvESSI/Connectivity",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(DataSpace);
		H5Dwrite(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Domain_Connectivity[domain_no]);
		H5Dclose(datasetId);

		dims1_out[0]= this->Domain_Pseudo_Number_of_Nodes[domain_no]; 
		DataSpace = H5Screate_simple(1, dims1_out, NULL);
		datasetId = H5Dcreate(id_File,"pvESSI/Inverse_Node_Map",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(DataSpace);
		H5Dwrite(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Domain_Inverse_Node_Map[domain_no]);
		H5Dclose(datasetId);

		dims1_out[0]= Domain_Pseudo_Number_of_Elements[domain_no]; 
		DataSpace = H5Screate_simple(1, dims1_out, NULL);
		datasetId = H5Dcreate(id_File,"pvESSI/Inverse_Element_Map",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(DataSpace);
		H5Dwrite(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Domain_Inverse_Element_Map[domain_no]);
		H5Dclose(datasetId);

		/************************************************* Write default Values of -1 for all time steps *****************************************/
		dims1_out[0]= this->Number_of_Time_Steps; 
		DataSpace = H5Screate_simple(1, dims1_out, NULL);
		datasetId = H5Dcreate(id_File,"pvESSI/Field_at_Nodes/Whether_Stress_Strain_Build",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(DataSpace);

		int *Whether_Stress_Strain_Build; Whether_Stress_Strain_Build = new int[Number_of_Time_Steps];
		int index = -1;
		while(++index<Number_of_Time_Steps){
			Whether_Stress_Strain_Build[index]=-1;
		}

		count1[0]   =Number_of_Time_Steps;
		dims1_out[0]=Number_of_Time_Steps;
		index_i     =0;

	    HDF5_Write_INT_Array_Data(datasetId, 1, dims1_out, &index_i, NULL, count1, NULL, Whether_Stress_Strain_Build); // Whether_Stress_Strain_Build
	    H5Dclose(datasetId);

		// *********************************************************** Creating Strain Dataset ****************************************************
		dims3[0]=this->Domain_Number_of_Nodes[domain_no];       dims3[1]= this->Number_of_Time_Steps;   dims3[2]= this->Number_of_Strain_Strain_Info; 
		maxdims3[0]=this->Domain_Number_of_Nodes[domain_no]; maxdims3[1]= this->Number_of_Time_Steps;   maxdims3[2]= this->Number_of_Strain_Strain_Info;
		DataSpace = H5Screate_simple(3, dims3, maxdims3);

	    hid_t prop = H5Pcreate (H5P_DATASET_CREATE);
	    hsize_t chunk_dims[3] = {this->Domain_Number_of_Nodes[domain_no],1,this->Number_of_Strain_Strain_Info};
	    H5Pset_chunk (prop, 3, chunk_dims);
		datasetId= H5Dcreate(id_File,"pvESSI/Field_at_Nodes/Stress_And_Strain",H5T_NATIVE_DOUBLE,DataSpace,H5P_DEFAULT,prop, H5P_DEFAULT); 
		status = H5Sclose(DataSpace);
		H5Dclose(datasetId);

		delete[] Whether_Stress_Strain_Build; Whether_Stress_Strain_Build = NULL;
	}
	else{
		cout << "<<<pvESSI>>>  File is not writable " << endl;
	}

	H5Fclose(this->id_File);
	this->id_File = H5Fopen(filename, id_H5F_READ_ONLY, id_H5F_CLOSE_STRONG);

}

//######################################################################
// Reads the Local Domain Maps if it has already been written in file. 
// inside pvESSI folder
void pvESSI::Read_Local_Domain_Maps(int domain_no){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Read_Local_Domain_Maps");
#endif

	if(!Domain_Basic_Info_Initialized[domain_no])
	{
		hid_t datasetId;

		// Finiding Pseudo number of nodes 
	   	datasetId = H5Dopen(id_File, "Model/Nodes/Number_of_DOFs", H5P_DEFAULT);
	    DataSpace = H5Dget_space(datasetId);
		H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
		this->Domain_Pseudo_Number_of_Nodes[domain_no] = dims1_out[0];
		H5Sclose(DataSpace);
		H5Dclose(datasetId);

		// Finiding Pseudo number of elements 
		datasetId = H5Dopen(id_File, "Model/Elements/Class_Tags", H5P_DEFAULT);
		DataSpace = H5Dget_space(datasetId);
		H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
		this->Domain_Pseudo_Number_of_Elements[domain_no] = dims1_out[0];
		H5Sclose(DataSpace);
		H5Dclose(datasetId);

		// Finiding Number of connectivity nodes
		datasetId = H5Dopen(id_File, "Model/Elements/Connectivity", H5P_DEFAULT);
		DataSpace = H5Dget_space(datasetId);
		H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
		this->Domain_Number_of_Connectivity_Nodes[domain_no] = dims1_out[0];
		H5Sclose(DataSpace);
		H5Dclose(datasetId);

		// Finding Number of Constrained nodes
		datasetId = H5Dopen(id_File, "/Model/Nodes/Constrained_DOFs", H5P_DEFAULT);
		DataSpace = H5Dget_space(datasetId);
		H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
		this->Domain_Number_of_Constrained_Dofs[domain_no] = dims1_out[0];
		H5Sclose(DataSpace);
		H5Dclose(datasetId);

		// Finiding total Number of gauss Points
		datasetId   = H5Dopen(id_File, "/Number_of_Gauss_Points", H5P_DEFAULT);
		H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Gauss_Points );
		this->Domain_Number_of_Gauss_Points[domain_no]=this->Number_of_Gauss_Points;
		H5Dclose(datasetId);

		// Finiding total Number of elements
		datasetId   = H5Dopen(id_File, "/Number_of_Elements", H5P_DEFAULT);
		H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Elements );
		this->Domain_Number_of_Elements[domain_no]=this->Number_of_Elements;
		H5Dclose(datasetId);

		// Finiding total Number of nodes
		datasetId  = H5Dopen(id_File, "/Number_of_Nodes", H5P_DEFAULT);
		H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Nodes );
		this->Domain_Number_of_Nodes[domain_no]=this->Number_of_Nodes;
		H5Dclose(datasetId);	

		datasetId                        = H5Dopen(id_File, "pvESSI/Node_Map", H5P_DEFAULT);
		INT_Time_Data_From_1_D_Dataset(datasetId, &Domain_Node_Map[domain_no]);
		H5Dclose(datasetId);

		datasetId                       = H5Dopen(id_File, "pvESSI/Element_Map", H5P_DEFAULT);
		INT_Time_Data_From_1_D_Dataset(datasetId, &Domain_Element_Map[domain_no]);
		H5Dclose(datasetId);

		datasetId                       = H5Dopen(id_File, "pvESSI/Number_of_DOfs", H5P_DEFAULT);
		INT_Time_Data_From_1_D_Dataset(datasetId, &Domain_Number_of_Dofs[domain_no]);
		H5Dclose(datasetId);

		datasetId                       = H5Dopen(id_File, "pvESSI/Number_of_Elements_Shared", H5P_DEFAULT);
		INT_Time_Data_From_1_D_Dataset(datasetId, &Domain_Number_of_Elements_Shared[domain_no]);
		H5Dclose(datasetId);

		datasetId                       = H5Dopen(id_File, "pvESSI/Number_of_Gauss_Elements_Shared", H5P_DEFAULT);
		INT_Time_Data_From_1_D_Dataset(datasetId, &Domain_Number_of_Gauss_Elements_Shared[domain_no]);
		H5Dclose(datasetId);

		datasetId                       = H5Dopen(id_File, "pvESSI/Class_Tags", H5P_DEFAULT);
		INT_Time_Data_From_1_D_Dataset(datasetId, &Domain_Class_Tags[domain_no]);
		H5Dclose(datasetId);

		datasetId                       = H5Dopen(id_File, "pvESSI/Connectivity", H5P_DEFAULT);
		INT_Time_Data_From_1_D_Dataset(datasetId, &Domain_Connectivity[domain_no]);
		H5Dclose(datasetId);

		datasetId                       = H5Dopen(id_File, "pvESSI/Inverse_Node_Map", H5P_DEFAULT);
		INT_Time_Data_From_1_D_Dataset(datasetId, &Domain_Inverse_Node_Map[domain_no]);
		H5Dclose(datasetId);

		datasetId                       = H5Dopen(id_File, "pvESSI/Inverse_Element_Map", H5P_DEFAULT);
		INT_Time_Data_From_1_D_Dataset(datasetId, &Domain_Inverse_Element_Map[domain_no]);
		H5Dclose(datasetId);

	}

	this->Domain_Basic_Info_Initialized[domain_no]=true;

}

/*****************************************************************************
* Builds time map i.e. a way to get the index number from a given non-interger
* time request from paraview desk. Although it works fine, its not the great 
* way to achieve it. 
* Needs to be changed. I think, paraview has fixed this issue in their latest 
* version, but I am not sure, need to check it. 
*****************************************************************************/
void pvESSI::Build_Time_Map(){


#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Time_Map");
#endif

	for (int i = 0;i<Number_of_Time_Steps ; i++){

		// check if already in map 
		// double store_time_step_no = 0;
		if(Time_Map.find(Time[i])== Time_Map.end()){
			Time_Map[Time[i]] = i;
			// store_time_step_no = i;
		}
		else{
			cout << "<<<<pvESSI>>>> Reappearance of time_step " << Time[i] << "s . Will show the visialization for original time" << endl;
		
		}
	}

	return;
}

/*****************************************************************************
* Updates the time vector information of the visualizaion
*****************************************************************************/
void pvESSI::Update_Time_Steps(){

	hid_t datasetId; 

	H5Eset_auto (NULL, NULL, NULL);  // To stop HDF5 from printing error message
	H5Fclose(this->id_File); // close if any instance of file is open
    /***************** File_id **********************************/
    this->id_File = H5Fopen(this->FileName, id_H5F_READ_ONLY, id_H5F_CLOSE_STRONG);;  

	/***************** Time Steps *******************************/
    datasetId = H5Dopen(id_File, "/Number_of_Time_Steps", H5P_DEFAULT);   
	H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Time_Steps);
	H5Dclose(datasetId);

	H5Eset_auto (NULL, NULL, NULL);  // To stop HDF5 from printing error message
	this->id_Eigen_Mode_Analysis = H5Gopen(id_File, "Eigen_Mode_Analysis", H5P_DEFAULT);
	if(this->id_Eigen_Mode_Analysis>0){
		datasetId  = H5Dopen(id_File,"Eigen_Mode_Analysis/Number_of_Eigen_Modes",H5P_DEFAULT);
		H5Dread(datasetId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Time_Steps);
		H5Dclose(datasetId);

		Number_of_Time_Steps = Number_of_Time_Steps +1;
		cout << "<<<<pvESSI>>>> Eigen_Mode_Analysis is On!!! \n" << endl;
		eigen_mode_on = true;
	}
	H5Gclose(this->id_Eigen_Mode_Analysis);

	float* Temp_Time;
	datasetId = H5Dopen(id_File, "/time", H5P_DEFAULT); 
	FLOAT_Time_Data_From_1_D_Dataset(datasetId,&Temp_Time);
	H5Dclose(datasetId);

	// close the file
	H5Fclose(this->id_File);

	Time = new double[this->Number_of_Time_Steps];

	if(eigen_mode_on)
	{
		// Initializing Time vector
		for(int p=0; p<Number_of_Time_Steps;p++)
			this->Time[p]=p;
	}
	else{
		// Initializing Time vector
		for(int p=0; p<Number_of_Time_Steps;p++)
			this->Time[p]=Temp_Time[p];
			// this->Time[p] = p;
	}

	delete [] Temp_Time; Temp_Time = NULL;
}

/*****************************************************************************
* Builds Meta Array map, so that it is easie to add attributes to vtkObject 
* and is neat and clean.
*****************************************************************************/
void pvESSI::Build_Meta_Array_Map(){


#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Meta_Array_Map");
#endif

	int key = 0;

	/*0  */ Meta_Array_Map["Generalized_Displacements"] = key; key=key+1;
	/*1  */ Meta_Array_Map["Generalized_Velocity"] = key; key=key+1;
	/*2  */ Meta_Array_Map["Generalized_Acceleration"] = key; key=key+1;
	/*3  */ Meta_Array_Map["Elastic_Strain"] = key; key=key+1;
	/*4  */ Meta_Array_Map["Plastic_Strain"] = key; key=key+1;
	/*5  */ Meta_Array_Map["Stress"] = key; key=key+1;
	/*6  */ Meta_Array_Map["Material_Tag"] = key; key=key+1;
	/*7  */ Meta_Array_Map["Total_Energy"] = key; key=key+1;
	/*8  */ Meta_Array_Map["Incremental_Energy"] = key; key=key+1;
	/*9  */ Meta_Array_Map["Node_Tag"] = key; key=key+1;
	/*10 */ Meta_Array_Map["Element_Tag"] = key; key=key+1;
	/*11 */ Meta_Array_Map["q"] = key; key=key+1;
	/*12 */ Meta_Array_Map["p"] = key; key=key+1;
	/*13 */ Meta_Array_Map["Plastic_Strain_q"] = key; key=key+1;
	/*14 */ Meta_Array_Map["Plastic_Strain_p"] = key; key=key+1;
	/*15 */ Meta_Array_Map["Support_Reactions"] = key; key=key+1;	
	/*16 */ Meta_Array_Map["Boundary_Conditions"] = key; key=key+1;
	/*17 */ Meta_Array_Map["Class_Tag"] = key; key=key+1;
	/*18 */ Meta_Array_Map["Partition_Info"] = key; key=key+1;
	/*19 */ Meta_Array_Map["Fluid_Displacements"] = key; key=key+1;
	/*20 */ Meta_Array_Map["Pore_Pressure"] = key; key=key+1;
}

void pvESSI::Build_Inverse_Matrices(){


#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Inverse_Matrices");
#endif

	Build_Brick_Coordinates();

	// for(int i=0; i<8; i++){
	// 	cout << Brick_8_Gauss_Coordinates[i][0] << " " << Brick_8_Gauss_Coordinates[i][1] << " " << Brick_8_Gauss_Coordinates[i][2] << endl;
	// }
	// cout << endl;

	double SQRT_1_3 = 1.0/sqrt(3.0);
	double SQRT_3_5 = sqrt(0.6);

	/////////////////////////////////////////// Building 8_Node_Brick_Inverse ///////////////////////////////////////

	this->Eight_Node_Brick_Inverse = new double*[8];

	double **Temp_Eight_Node_Brick = new double*[8];
	for(int i = 0; i < 8; ++i){
    	Temp_Eight_Node_Brick[i] = new double[8];
    	Eight_Node_Brick_Inverse[i] = new double[8];
	}

	for (int j=0; j<8;j++){ // Gauss Point Coordinates
		for (int i=0; i<8;i++){ // Shape Function
			Temp_Eight_Node_Brick[j][i] =   0.125*(1+Brick_Coordinates[i][0]*Brick_8_Gauss_Coordinates[j][0]*SQRT_1_3)*
											      (1+Brick_Coordinates[i][1]*Brick_8_Gauss_Coordinates[j][1]*SQRT_1_3)*
											      (1+Brick_Coordinates[i][2]*Brick_8_Gauss_Coordinates[j][2]*SQRT_1_3);
		}
	}

	vtkMath::InvertMatrix(Temp_Eight_Node_Brick,Eight_Node_Brick_Inverse,8);

	////////////////////////////////////// Building 20_Node_Brick_Inverse //////////////////////////////////////

  	this->Twenty_Node_Brick_Inverse = new double*[20];
	double **Temp_Twenty_Node_Brick = new double*[20];
	for(int i = 0; i < 20; ++i){
    	Temp_Twenty_Node_Brick[i] = new double[20];
    	Twenty_Node_Brick_Inverse[i] = new double[20];
	}
	for (int j=0; j<20;j++){ // Gauss Point Coordinates 
		for (int i=0; i<8;i++){ // Shape Function 
			Temp_Twenty_Node_Brick[j][i] =   0.125*(1+Brick_Coordinates[i][0]*Brick_20_Gauss_Coordinates[j][0]*SQRT_3_5)* 
												   (1+Brick_Coordinates[i][1]*Brick_20_Gauss_Coordinates[j][1]*SQRT_3_5)* 
												   (1+Brick_Coordinates[i][2]*Brick_20_Gauss_Coordinates[j][2]*SQRT_3_5)* 
												   (-2+ (Brick_Coordinates[i][0]*Brick_20_Gauss_Coordinates[j][0]*SQRT_3_5)+ 
												        (Brick_Coordinates[i][1]*Brick_20_Gauss_Coordinates[j][1]*SQRT_3_5)+ 
												        (Brick_Coordinates[i][2]*Brick_20_Gauss_Coordinates[j][2]*SQRT_3_5));
		}
		int x[] = { 8,10,12,14};
		for (int k=0; k<4;k++){ // Shape Function 
			int i = x[k];
			Temp_Twenty_Node_Brick[j][i] = 	  0.25*(1-pow(Brick_20_Gauss_Coordinates[j][0]*SQRT_3_5,2))*
												   (1+Brick_Coordinates[i][1]*Brick_20_Gauss_Coordinates[j][1]*SQRT_3_5)*
												   (1+Brick_Coordinates[i][2]*Brick_20_Gauss_Coordinates[j][2]*SQRT_3_5);
		}
		x[0] = 9; x[1]=11; x[2]=13; x[3]=15;
		for (int k=0; k<4;k++){ // Shape Function 
			int i = x[k];
			Temp_Twenty_Node_Brick[j][i] =    0.25*(1-pow(Brick_20_Gauss_Coordinates[j][1]*SQRT_3_5,2))*
												   (1+Brick_Coordinates[i][0]*Brick_20_Gauss_Coordinates[j][0]*SQRT_3_5)*
												   (1+Brick_Coordinates[i][2]*Brick_20_Gauss_Coordinates[j][2]*SQRT_3_5);
		}
		x[0]=16; x[1]=17; x[2]=18; x[3]=19;
		for (int k=0; k<4;k++){ // Shape Function 
			int i = x[k];
			Temp_Twenty_Node_Brick[j][i] =    0.25*(1-pow(Brick_20_Gauss_Coordinates[j][2]*SQRT_3_5,2))*
												   (1+Brick_Coordinates[i][0]*Brick_20_Gauss_Coordinates[j][0]*SQRT_3_5)*
												   (1+Brick_Coordinates[i][1]*Brick_20_Gauss_Coordinates[j][1]*SQRT_3_5);
		}

	}

	vtkMath::InvertMatrix(Temp_Twenty_Node_Brick,Twenty_Node_Brick_Inverse,20);

	// for (int j=0; j<20;j++){
	// 	for (int i=0; i<20;i++){
	// 		cout << std::setw(8) << std::setprecision(5) << Twenty_Node_Brick_Inverse[j][i] << " " ;
	// 	}
	// 	cout << endl;
	// }

	////////////////////////////////////////// Building 27_Node_Brick_Inverse //////////////////////////////////////

  	this->Twenty_Seven_Node_Brick_Inverse = new double*[27];
	double **Temp_Twenty_Seven_Brick = new double*[27];
	for(int i = 0; i < 27; ++i){
    	Temp_Twenty_Seven_Brick[i] = new double[27];
    	Twenty_Seven_Node_Brick_Inverse[i] = new double[27];
	}
	for (int j=0; j<27;j++){
		for (int i=0; i<8;i++){
			Temp_Twenty_Seven_Brick[j][i]         = 0.125*(1+Brick_Coordinates[i][0]*Brick_27_Gauss_Coordinates[j][0]*SQRT_3_5)*
												          (1+Brick_Coordinates[i][1]*Brick_27_Gauss_Coordinates[j][1]*SQRT_3_5)*
												          (1+Brick_Coordinates[i][2]*Brick_27_Gauss_Coordinates[j][2]*SQRT_3_5)*
												          ( (Brick_Coordinates[i][0]*Brick_27_Gauss_Coordinates[j][0]*SQRT_3_5)*
												            (Brick_Coordinates[i][1]*Brick_27_Gauss_Coordinates[j][1]*SQRT_3_5)*
												            (Brick_Coordinates[i][2]*Brick_27_Gauss_Coordinates[j][2]*SQRT_3_5));
		}
		int x[] = { 8,10,12,14};
		for (int k=0; k<4;k++){
			int i = x[k];
			Temp_Twenty_Seven_Brick[j][i]         = 0.25*(1-pow(Brick_27_Gauss_Coordinates[j][0]*SQRT_3_5,2))*
												         (1+Brick_Coordinates[i][1]*Brick_27_Gauss_Coordinates[j][1]*SQRT_3_5)*
												         (1+Brick_Coordinates[i][2]*Brick_27_Gauss_Coordinates[j][2]*SQRT_3_5)*
												         ( (Brick_Coordinates[i][1]*Brick_27_Gauss_Coordinates[j][1]*SQRT_3_5)*
												           (Brick_Coordinates[i][2]*Brick_27_Gauss_Coordinates[j][2]*SQRT_3_5));
		}
		x[0] = 9; x[1]=11; x[2]=13; x[3]=15;
		for (int k=0; k<4;k++){
			int i = x[k];
			Temp_Twenty_Seven_Brick[j][i]         = 0.25*(1-pow(Brick_27_Gauss_Coordinates[j][1]*SQRT_3_5,2))*
												         (1+Brick_Coordinates[i][0]*Brick_27_Gauss_Coordinates[j][0]*SQRT_3_5)*
												         (1+Brick_Coordinates[i][2]*Brick_27_Gauss_Coordinates[j][2]*SQRT_3_5)*
												         ( (Brick_Coordinates[i][0]*Brick_27_Gauss_Coordinates[j][0]*SQRT_3_5)*
												           (Brick_Coordinates[i][2]*Brick_27_Gauss_Coordinates[j][2]*SQRT_3_5));
		}
		x[0]=16; x[1]=17; x[2]=18; x[3]=19;
		for (int k=0; k<4;k++){
			int i = x[k];
			Temp_Twenty_Seven_Brick[j][i]         = 0.25*(1-pow(Brick_27_Gauss_Coordinates[j][2]*SQRT_3_5,2))*
												         (1+Brick_Coordinates[i][0]*Brick_27_Gauss_Coordinates[j][0]*SQRT_3_5)*
												         (1+Brick_Coordinates[i][1]*Brick_27_Gauss_Coordinates[j][1]*SQRT_3_5)*
												         ( (Brick_Coordinates[i][0]*Brick_27_Gauss_Coordinates[j][0]*SQRT_3_5)*
												           (Brick_Coordinates[i][1]*Brick_27_Gauss_Coordinates[j][1]*SQRT_3_5));
		}
		{
			int i =20;
			Temp_Twenty_Seven_Brick[j][i]         = 0.25*(1-pow(Brick_27_Gauss_Coordinates[j][0]*SQRT_3_5,2))*
													     (1-pow(Brick_27_Gauss_Coordinates[j][1]*SQRT_3_5,2))*
													     (1-pow(Brick_27_Gauss_Coordinates[j][2]*SQRT_3_5,2));
		}
		int y[] = {21,23};
		for (int k=0; k<2;k++){
			int i = y[k];
			Temp_Twenty_Seven_Brick[j][i]         = 0.5*(1-pow(Brick_27_Gauss_Coordinates[j][0]*SQRT_3_5,2))*
													    (1-pow(Brick_27_Gauss_Coordinates[j][2]*SQRT_3_5,2))*
													    (1+Brick_Coordinates[i][1]*Brick_27_Gauss_Coordinates[j][1]*SQRT_3_5)*
													      (Brick_Coordinates[i][1]*Brick_27_Gauss_Coordinates[j][1]*SQRT_3_5);
		}
		y[0]=22; y[1]=24;
		for (int k=0; k<2;k++){
			int i = y[k];
			Temp_Twenty_Seven_Brick[j][i]         = 0.5*(1-pow(Brick_27_Gauss_Coordinates[j][1]*SQRT_3_5,2))*
													    (1-pow(Brick_27_Gauss_Coordinates[j][2]*SQRT_3_5,2))*
													    (1+Brick_Coordinates[i][0]*Brick_27_Gauss_Coordinates[j][0]*SQRT_3_5)*
													      (Brick_Coordinates[i][0]*Brick_27_Gauss_Coordinates[j][0]*SQRT_3_5);
		}
		y[0]=25; y[1]=26;
		for (int k=0; k<2;k++){
			int i = y[k];
			Temp_Twenty_Seven_Brick[j][i]         = 0.5*(1-pow(Brick_27_Gauss_Coordinates[j][0]*SQRT_3_5,2))*
													    (1-pow(Brick_27_Gauss_Coordinates[j][1]*SQRT_3_5,2))*
													    (1+Brick_Coordinates[i][2]*Brick_27_Gauss_Coordinates[j][2]*SQRT_3_5)*
													      (Brick_Coordinates[i][2]*Brick_27_Gauss_Coordinates[j][2]*SQRT_3_5);
		}
	}

	vtkMath::InvertMatrix(Temp_Twenty_Seven_Brick,Twenty_Seven_Node_Brick_Inverse,27);

	// ///////////////////////////// Printing For Debugging //////////////////////////////////////////////////////////
	// cout << " ********************************************************************* " << endl << endl;
	// for (int j=0; j<27;j++){
	// 	for (int i=0; i<27;i++){
	// 		cout << std::setw(8) << std::setprecision(5) << Twenty_Seven_Node_Brick_Inverse[j][i] << " " ;
	// 	}

	// 	cout << endl;
	// }

	delete [] Temp_Twenty_Seven_Brick;Temp_Twenty_Seven_Brick=NULL;
	delete [] Temp_Twenty_Node_Brick;Temp_Twenty_Node_Brick=NULL;
	delete [] Temp_Eight_Node_Brick;Temp_Eight_Node_Brick=NULL;

	return;
  }

void pvESSI::Build_Brick_Coordinates(){


#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Brick_Coordinates");
#endif


	Brick_8_Gauss_Coordinates[0 ][0]=-1;  Brick_8_Gauss_Coordinates[0 ][1]=-1; Brick_8_Gauss_Coordinates[0 ][2]=-1;
	Brick_8_Gauss_Coordinates[1 ][0]=-1;  Brick_8_Gauss_Coordinates[1 ][1]=-1; Brick_8_Gauss_Coordinates[1 ][2]= 1;
	Brick_8_Gauss_Coordinates[2 ][0]=-1;  Brick_8_Gauss_Coordinates[2 ][1]= 1; Brick_8_Gauss_Coordinates[2 ][2]=-1;
	Brick_8_Gauss_Coordinates[3 ][0]=-1;  Brick_8_Gauss_Coordinates[3 ][1]= 1; Brick_8_Gauss_Coordinates[3 ][2]= 1;
	Brick_8_Gauss_Coordinates[4 ][0]= 1;  Brick_8_Gauss_Coordinates[4 ][1]=-1; Brick_8_Gauss_Coordinates[4 ][2]=-1;
	Brick_8_Gauss_Coordinates[5 ][0]= 1;  Brick_8_Gauss_Coordinates[5 ][1]=-1; Brick_8_Gauss_Coordinates[5 ][2]= 1;
	Brick_8_Gauss_Coordinates[6 ][0]= 1;  Brick_8_Gauss_Coordinates[6 ][1]= 1; Brick_8_Gauss_Coordinates[6 ][2]=-1;
	Brick_8_Gauss_Coordinates[7 ][0]= 1;  Brick_8_Gauss_Coordinates[7 ][1]= 1; Brick_8_Gauss_Coordinates[7 ][2]= 1;


	Brick_20_Gauss_Coordinates[0 ][0]=-1;  Brick_20_Gauss_Coordinates[0 ][1]=-1; Brick_20_Gauss_Coordinates[0 ][2]=-1;
	Brick_20_Gauss_Coordinates[1 ][0]=-1;  Brick_20_Gauss_Coordinates[1 ][1]=-1; Brick_20_Gauss_Coordinates[1 ][2]= 0;
	Brick_20_Gauss_Coordinates[2 ][0]=-1;  Brick_20_Gauss_Coordinates[2 ][1]=-1; Brick_20_Gauss_Coordinates[2 ][2]= 1;
	Brick_20_Gauss_Coordinates[3 ][0]=-1;  Brick_20_Gauss_Coordinates[3 ][1]= 0; Brick_20_Gauss_Coordinates[3 ][2]=-1;
	Brick_20_Gauss_Coordinates[4 ][0]=-1;  Brick_20_Gauss_Coordinates[4 ][1]= 0; Brick_20_Gauss_Coordinates[4 ][2]= 1;
	Brick_20_Gauss_Coordinates[5 ][0]=-1;  Brick_20_Gauss_Coordinates[5 ][1]= 1; Brick_20_Gauss_Coordinates[5 ][2]=-1;
	Brick_20_Gauss_Coordinates[6 ][0]=-1;  Brick_20_Gauss_Coordinates[6 ][1]= 1; Brick_20_Gauss_Coordinates[6 ][2]= 0;
	Brick_20_Gauss_Coordinates[7 ][0]=-1;  Brick_20_Gauss_Coordinates[7 ][1]= 1; Brick_20_Gauss_Coordinates[7 ][2]= 1;
	Brick_20_Gauss_Coordinates[8 ][0]= 0;  Brick_20_Gauss_Coordinates[8 ][1]=-1; Brick_20_Gauss_Coordinates[8 ][2]=-1;
	Brick_20_Gauss_Coordinates[9 ][0]= 0;  Brick_20_Gauss_Coordinates[9 ][1]=-1; Brick_20_Gauss_Coordinates[9 ][2]= 1;
	Brick_20_Gauss_Coordinates[10][0]= 0;  Brick_20_Gauss_Coordinates[10][1]= 1; Brick_20_Gauss_Coordinates[10][2]=-1;
	Brick_20_Gauss_Coordinates[11][0]= 0;  Brick_20_Gauss_Coordinates[11][1]= 1; Brick_20_Gauss_Coordinates[11][2]= 1;
	Brick_20_Gauss_Coordinates[12][0]= 1;  Brick_20_Gauss_Coordinates[12][1]=-1; Brick_20_Gauss_Coordinates[12][2]=-1;
	Brick_20_Gauss_Coordinates[13][0]= 1;  Brick_20_Gauss_Coordinates[13][1]=-1; Brick_20_Gauss_Coordinates[13][2]= 0;
	Brick_20_Gauss_Coordinates[14][0]= 1;  Brick_20_Gauss_Coordinates[14][1]=-1; Brick_20_Gauss_Coordinates[14][2]= 1;
	Brick_20_Gauss_Coordinates[15][0]= 1;  Brick_20_Gauss_Coordinates[15][1]= 0; Brick_20_Gauss_Coordinates[15][2]=-1;
	Brick_20_Gauss_Coordinates[16][0]= 1;  Brick_20_Gauss_Coordinates[16][1]= 0; Brick_20_Gauss_Coordinates[16][2]= 1;
	Brick_20_Gauss_Coordinates[17][0]= 1;  Brick_20_Gauss_Coordinates[17][1]= 1; Brick_20_Gauss_Coordinates[17][2]=-1;
	Brick_20_Gauss_Coordinates[18][0]= 1;  Brick_20_Gauss_Coordinates[18][1]= 1; Brick_20_Gauss_Coordinates[18][2]= 0;
	Brick_20_Gauss_Coordinates[19][0]= 1;  Brick_20_Gauss_Coordinates[19][1]= 1; Brick_20_Gauss_Coordinates[19][2]= 1;

	Brick_27_Gauss_Coordinates[0 ][0]=-1;  Brick_27_Gauss_Coordinates[0 ][1]=-1; Brick_27_Gauss_Coordinates[0 ][2]=-1;
	Brick_27_Gauss_Coordinates[1 ][0]=-1;  Brick_27_Gauss_Coordinates[1 ][1]=-1; Brick_27_Gauss_Coordinates[1 ][2]= 0;
	Brick_27_Gauss_Coordinates[2 ][0]=-1;  Brick_27_Gauss_Coordinates[2 ][1]=-1; Brick_27_Gauss_Coordinates[2 ][2]= 1;
	Brick_27_Gauss_Coordinates[3 ][0]=-1;  Brick_27_Gauss_Coordinates[3 ][1]= 0; Brick_27_Gauss_Coordinates[3 ][2]=-1;
	Brick_27_Gauss_Coordinates[4 ][0]=-1;  Brick_27_Gauss_Coordinates[4 ][1]= 0; Brick_27_Gauss_Coordinates[4 ][2]= 0;
	Brick_27_Gauss_Coordinates[5 ][0]=-1;  Brick_27_Gauss_Coordinates[5 ][1]= 0; Brick_27_Gauss_Coordinates[5 ][2]= 1;
	Brick_27_Gauss_Coordinates[6 ][0]=-1;  Brick_27_Gauss_Coordinates[6 ][1]= 1; Brick_27_Gauss_Coordinates[6 ][2]=-1;
	Brick_27_Gauss_Coordinates[7 ][0]=-1;  Brick_27_Gauss_Coordinates[7 ][1]= 1; Brick_27_Gauss_Coordinates[7 ][2]= 0;
	Brick_27_Gauss_Coordinates[8 ][0]=-1;  Brick_27_Gauss_Coordinates[8 ][1]= 1; Brick_27_Gauss_Coordinates[8 ][2]= 1;
	Brick_27_Gauss_Coordinates[9 ][0]= 0;  Brick_27_Gauss_Coordinates[9 ][1]=-1; Brick_27_Gauss_Coordinates[9 ][2]=-1;
	Brick_27_Gauss_Coordinates[10][0]= 0;  Brick_27_Gauss_Coordinates[10][1]=-1; Brick_27_Gauss_Coordinates[10][2]= 0;
	Brick_27_Gauss_Coordinates[11][0]= 0;  Brick_27_Gauss_Coordinates[11][1]=-1; Brick_27_Gauss_Coordinates[11][2]= 1;
	Brick_27_Gauss_Coordinates[12][0]= 0;  Brick_27_Gauss_Coordinates[12][1]= 0; Brick_27_Gauss_Coordinates[12][2]=-1;
	Brick_27_Gauss_Coordinates[13][0]= 0;  Brick_27_Gauss_Coordinates[13][1]= 0; Brick_27_Gauss_Coordinates[13][2]= 0;
	Brick_27_Gauss_Coordinates[14][0]= 0;  Brick_27_Gauss_Coordinates[14][1]= 0; Brick_27_Gauss_Coordinates[14][2]= 1;
	Brick_27_Gauss_Coordinates[15][0]= 0;  Brick_27_Gauss_Coordinates[15][1]= 1; Brick_27_Gauss_Coordinates[15][2]=-1;
	Brick_27_Gauss_Coordinates[16][0]= 0;  Brick_27_Gauss_Coordinates[16][1]= 1; Brick_27_Gauss_Coordinates[16][2]= 0;
	Brick_27_Gauss_Coordinates[17][0]= 0;  Brick_27_Gauss_Coordinates[17][1]= 1; Brick_27_Gauss_Coordinates[17][2]= 1;
	Brick_27_Gauss_Coordinates[18][0]= 1;  Brick_27_Gauss_Coordinates[18][1]=-1; Brick_27_Gauss_Coordinates[18][2]=-1;
	Brick_27_Gauss_Coordinates[19][0]= 1;  Brick_27_Gauss_Coordinates[19][1]=-1; Brick_27_Gauss_Coordinates[19][2]= 0;
	Brick_27_Gauss_Coordinates[20][0]= 1;  Brick_27_Gauss_Coordinates[20][1]=-1; Brick_27_Gauss_Coordinates[20][2]= 1;
	Brick_27_Gauss_Coordinates[21][0]= 1;  Brick_27_Gauss_Coordinates[21][1]= 0; Brick_27_Gauss_Coordinates[21][2]=-1;
	Brick_27_Gauss_Coordinates[22][0]= 1;  Brick_27_Gauss_Coordinates[22][1]= 0; Brick_27_Gauss_Coordinates[22][2]= 0;
	Brick_27_Gauss_Coordinates[23][0]= 1;  Brick_27_Gauss_Coordinates[23][1]= 0; Brick_27_Gauss_Coordinates[23][2]= 1;
	Brick_27_Gauss_Coordinates[24][0]= 1;  Brick_27_Gauss_Coordinates[24][1]= 1; Brick_27_Gauss_Coordinates[24][2]=-1;
	Brick_27_Gauss_Coordinates[25][0]= 1;  Brick_27_Gauss_Coordinates[25][1]= 1; Brick_27_Gauss_Coordinates[25][2]= 0;
	Brick_27_Gauss_Coordinates[26][0]= 1;  Brick_27_Gauss_Coordinates[26][1]= 1; Brick_27_Gauss_Coordinates[26][2]= 1;

	Brick_Coordinates[0 ][0]= 1;  Brick_Coordinates[0 ][1]= 1; Brick_Coordinates[0 ][2]= 1;
	Brick_Coordinates[1 ][0]=-1;  Brick_Coordinates[1 ][1]= 1; Brick_Coordinates[1 ][2]= 1;
	Brick_Coordinates[2 ][0]=-1;  Brick_Coordinates[2 ][1]=-1; Brick_Coordinates[2 ][2]= 1;
	Brick_Coordinates[3 ][0]= 1;  Brick_Coordinates[3 ][1]=-1; Brick_Coordinates[3 ][2]= 1;
	Brick_Coordinates[4 ][0]= 1;  Brick_Coordinates[4 ][1]= 1; Brick_Coordinates[4 ][2]=-1;
	Brick_Coordinates[5 ][0]=-1;  Brick_Coordinates[5 ][1]= 1; Brick_Coordinates[5 ][2]=-1;
	Brick_Coordinates[6 ][0]=-1;  Brick_Coordinates[6 ][1]=-1; Brick_Coordinates[6 ][2]=-1;
	Brick_Coordinates[7 ][0]= 1;  Brick_Coordinates[7 ][1]=-1; Brick_Coordinates[7 ][2]=-1;
	Brick_Coordinates[8 ][0]= 0;  Brick_Coordinates[8 ][1]= 1; Brick_Coordinates[8 ][2]= 1;
	Brick_Coordinates[9 ][0]=-1;  Brick_Coordinates[9 ][1]= 0; Brick_Coordinates[9 ][2]= 1;
	Brick_Coordinates[10][0]= 0;  Brick_Coordinates[10][1]=-1; Brick_Coordinates[10][2]= 1;
	Brick_Coordinates[11][0]= 1;  Brick_Coordinates[11][1]= 0; Brick_Coordinates[11][2]= 1;
	Brick_Coordinates[12][0]= 0;  Brick_Coordinates[12][1]= 1; Brick_Coordinates[12][2]=-1;
	Brick_Coordinates[13][0]=-1;  Brick_Coordinates[13][1]= 0; Brick_Coordinates[13][2]=-1;
	Brick_Coordinates[14][0]= 0;  Brick_Coordinates[14][1]=-1; Brick_Coordinates[14][2]=-1;
	Brick_Coordinates[15][0]= 1;  Brick_Coordinates[15][1]= 0; Brick_Coordinates[15][2]=-1;
	Brick_Coordinates[16][0]= 1;  Brick_Coordinates[16][1]= 1; Brick_Coordinates[16][2]= 0;
	Brick_Coordinates[17][0]=-1;  Brick_Coordinates[17][1]= 1; Brick_Coordinates[17][2]= 0;
	Brick_Coordinates[18][0]=-1;  Brick_Coordinates[18][1]=-1; Brick_Coordinates[18][2]= 0;
	Brick_Coordinates[19][0]= 1;  Brick_Coordinates[19][1]=-1; Brick_Coordinates[19][2]= 0;
	Brick_Coordinates[20][0]= 0;  Brick_Coordinates[20][1]= 0; Brick_Coordinates[20][2]= 0;
	Brick_Coordinates[21][0]= 0;  Brick_Coordinates[21][1]= 1; Brick_Coordinates[21][2]= 0;
	Brick_Coordinates[22][0]=-1;  Brick_Coordinates[22][1]= 0; Brick_Coordinates[22][2]= 0;
	Brick_Coordinates[23][0]= 0;  Brick_Coordinates[23][1]=-1; Brick_Coordinates[23][2]= 0;
	Brick_Coordinates[24][0]= 1;  Brick_Coordinates[24][1]= 0; Brick_Coordinates[24][2]= 0;
	Brick_Coordinates[25][0]= 0;  Brick_Coordinates[25][1]= 0; Brick_Coordinates[25][2]= 1;
	Brick_Coordinates[26][0]= 0;  Brick_Coordinates[26][1]= 0; Brick_Coordinates[26][2]=-1;

	return;

}


void pvESSI::Build_Gauss_To_Node_Interpolation_Map(){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Gauss_To_Node_Interpolation_Map");
#endif

	Gauss_To_Node_Interpolation_Map[2]  = this->Eight_Node_Brick_Inverse;
	Gauss_To_Node_Interpolation_Map[3]  = this->Eight_Node_Brick_Inverse;
	Gauss_To_Node_Interpolation_Map[4]  = this->Eight_Node_Brick_Inverse;
	Gauss_To_Node_Interpolation_Map[5]  = this->Twenty_Node_Brick_Inverse;
	Gauss_To_Node_Interpolation_Map[6]  = this->Twenty_Node_Brick_Inverse;
	Gauss_To_Node_Interpolation_Map[7]  = this->Twenty_Node_Brick_Inverse;
	Gauss_To_Node_Interpolation_Map[8]  = this->Twenty_Seven_Node_Brick_Inverse;
	Gauss_To_Node_Interpolation_Map[9]  = this->Twenty_Seven_Node_Brick_Inverse;
	Gauss_To_Node_Interpolation_Map[10] = this->Twenty_Seven_Node_Brick_Inverse;
}

void pvESSI::HDF5_Read_INT_Array_Data(hid_t id_DataSet, int rank, hsize_t *data_dims, hsize_t *offset, hsize_t *stride, hsize_t *count, hsize_t *block, int* data){

  //Get pointer to the dataspace and create the memory space
  DataSpace = H5Dget_space(id_DataSet);
  MemSpace  = H5Screate_simple(1, data_dims, data_dims);

  //Select the region of data to output to
  H5Sselect_hyperslab(
      DataSpace,             // Id of the parent dataspace
      H5S_SELECT_SET,        // Selection operatior H5S_SELECT_<>, where <> = {SET, OR, AND, XOR, NOTB, NOTA}
      offset,                // start of selection
      stride,                // stride in each dimension, NULL  is select everything
      count ,                // how many blocks to select in each direction
      block                  // little block selected per selection
  );

  H5Dread(
	      id_DataSet,              // Dataset to write to
	      H5T_NATIVE_INT,     // Format of data in memory
	      MemSpace,           // Description of data in memory
	      DataSpace,          // Description of data in storage (including selection)
	      H5P_DEFAULT,         // Form of reading
	      data                // The actual data
	  );

  H5Sclose(DataSpace);
  H5Sclose(MemSpace);
}


void pvESSI::HDF5_Write_INT_Array_Data(hid_t id_DataSet, int rank, hsize_t *data_dims, hsize_t *offset, hsize_t *stride, hsize_t *count, hsize_t *block, int* data){

  //Get pointer to the dataspace and create the memory space
  DataSpace = H5Dget_space(id_DataSet);
  MemSpace  = H5Screate_simple(1, data_dims, data_dims);

  //Select the region of data to output to
  H5Sselect_hyperslab(
      DataSpace,             // Id of the parent dataspace
      H5S_SELECT_SET,        // Selection operatior H5S_SELECT_<>, where <> = {SET, OR, AND, XOR, NOTB, NOTA}
      offset,                // start of selection
      stride,                // stride in each dimension, NULL  is select everything
      count ,                // how many blocks to select in each direction
      block                  // little block selected per selection
  );

  H5Dwrite(
	      id_DataSet,              // Dataset to write to
	      H5T_NATIVE_INT,     // Format of data in memory
	      MemSpace,           // Description of data in memory
	      DataSpace,          // Description of data in storage (including selection)
	      H5P_DEFAULT,         // Form of reading
	      data                // The actual data
	  );

  H5Sclose(DataSpace);
  H5Sclose(MemSpace);
}

void pvESSI::INTERPOLATE_FLOAT_Time_Data_From_2_D_Dataset(hid_t datasetId, int timeIndex1, int timeIndex2, float shapeFunction1, float shapeFunction2, float** DataArray){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::INTERPOLATE_FLOAT_Time_Data_From_2_D_Dataset");
	cout << "DatasetId            " << datasetId << endl;
#endif

	if(timeIndex1==timeIndex2 or shapeFunction1==1){
		 float *DataArray1; FLOAT_Time_Data_From_2_D_Dataset(datasetId, timeIndex1, &DataArray1);
		 *DataArray = DataArray1;
		 return;
	}
	else if(shapeFunction2==1){
		float *DataArray2; FLOAT_Time_Data_From_2_D_Dataset(datasetId, timeIndex2, &DataArray2);
		*DataArray = DataArray2;
		return;
	}

	float *DataArray1; FLOAT_Time_Data_From_2_D_Dataset(datasetId, timeIndex1, &DataArray1);
	float *DataArray2; FLOAT_Time_Data_From_2_D_Dataset(datasetId, timeIndex2, &DataArray2);

	*DataArray = new float[dims2_out[0]];
	for (int i=0;i<dims2_out[0];i++){
		(*DataArray)[i]=shapeFunction1*DataArray1[i]+shapeFunction2*DataArray2[i];
	}

	delete []  DataArray1; DataArray1 = NULL;
	delete []  DataArray2; DataArray2 = NULL;

}


void pvESSI::INTERPOLATE_FLOAT_Time_Data_From_3_D_Dataset(hid_t datasetId, int timeIndex1, int timeIndex2, float shapeFunction1, float shapeFunction2, float** DataArray){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::INTERPOLATE_FLOAT_Time_Data_From_3_D_Dataset");
	cout << "DatasetId            " << datasetId << endl;
#endif

	H5Drefresh(datasetId);
	DataSpace = H5Dget_space(datasetId);
	H5Sget_simple_extent_dims(DataSpace, dims3_out, NULL);
	float *DataArray1 = new float[dims3_out[0]*dims3_out[2]];
	offset3[0] = 0;  				 offset3[1] =timeIndex1;          offset3[2] = 0;
    count3 [0] = dims3_out[0];		 count3 [1] = 1;		    	  count3 [2] = dims3_out[2];
    dims2_out[0] = dims3_out[0];;	     dims2_out[1] = dims3_out[2];
    MemSpace = H5Screate_simple(2,dims2_out,NULL);
    H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset3,NULL,count3,NULL);
    H5Dread(datasetId, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, DataArray1); 
    H5Sclose(MemSpace);
    H5Sclose(DataSpace);

    H5Drefresh(datasetId);
	DataSpace = H5Dget_space(datasetId);
	H5Sget_simple_extent_dims(DataSpace, dims3_out, NULL);
	float *DataArray2 = new float[dims3_out[0]*dims3_out[2]];
	offset3[0] = 0;  				 offset3[1] =timeIndex1;          offset3[2] = 0;
    count3 [0] = dims3_out[0];		 count3 [1] = 1;		    	  count3 [2] = dims3_out[2];
    dims2_out[0] = dims3_out[0];;	     dims2_out[1] = dims3_out[2];
    MemSpace = H5Screate_simple(2,dims2_out,NULL);
    H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset3,NULL,count3,NULL);
    H5Dread(datasetId, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, DataArray2); 
    H5Sclose(MemSpace);
    H5Sclose(DataSpace);

	*DataArray = new float[dims3_out[0]*dims3_out[2]];
	for (int i=0;i<dims3_out[0]*dims3_out[2];i++){
		(*DataArray)[i]=shapeFunction1*DataArray1[i]+shapeFunction2*DataArray2[i];
	}

	delete []  DataArray1; DataArray1 = NULL;
	delete []  DataArray2; DataArray2 = NULL;

}


void pvESSI::FLOAT_Time_Data_From_2_D_Dataset(hid_t datasetId, int timeIndex, float** DataArray){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::FLOAT_Time_Data_From_2_D_Dataset");
	cout << "DatasetId            " << datasetId << endl;
#endif

	H5Drefresh(datasetId);
	DataSpace = H5Dget_space(datasetId);
	H5Sget_simple_extent_dims(DataSpace, dims2_out, NULL);	
	*DataArray = new float[dims2_out[0]];	
	offset2[0]=0; 			      count2[0] = dims2_out[0];		dims1_out[0]=dims2_out[0];
	offset2[1]=timeIndex; 		  count2[1] = 1;				MemSpace = H5Screate_simple(1,dims1_out,NULL);
	H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset2,NULL,count2,NULL);
	H5Dread(datasetId, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, *DataArray); 
	H5Sclose(MemSpace);
	H5Sclose(DataSpace); 

}


void pvESSI::FLOAT_Time_Data_From_3_D_Dataset(hid_t datasetId, int timeIndex, float** DataArray){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::FLOAT_Time_Data_From_3_D_Dataset");
	cout << "DatasetId            " << datasetId << endl;
#endif

	H5Drefresh(datasetId);
	DataSpace = H5Dget_space(datasetId);
	H5Sget_simple_extent_dims(DataSpace, dims3_out, NULL);
	*DataArray = new float[dims3_out[0]*dims3_out[2]];
	offset3[0] = 0;  				 offset3[1] =timeIndex;          offset3[2] = 0;
    count3 [0] = dims3_out[0];		 count3 [1] = 1;		    	  count3 [2] = dims3_out[2];
    dims2_out[0] = dims3_out[0];;	     dims2_out[1] = dims3_out[2];
    MemSpace = H5Screate_simple(2,dims2_out,NULL);
    H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset3,NULL,count3,NULL);
    H5Dread(datasetId, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT,*DataArray); 
    H5Sclose(MemSpace);
    H5Sclose(DataSpace);

}


void pvESSI::INT_Time_Data_From_1_D_Dataset(hid_t datasetId, int** DataArray){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::INT_Time_Data_From_1_D_Dataset");
	cout << "DatasetId            " << datasetId << endl;
#endif

	H5Drefresh(datasetId);
	DataSpace = H5Dget_space(datasetId);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
	*DataArray= new int[dims1_out[0]];
	offset1[0]=0;   	MemSpace = H5Screate_simple(1,dims1_out,NULL);
	H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset1,NULL,dims1_out,NULL);
	H5Dread(datasetId, H5T_NATIVE_INT, MemSpace, DataSpace, H5P_DEFAULT, *DataArray); 
	H5Sclose(MemSpace); status=H5Sclose(DataSpace);

}


void pvESSI::FLOAT_Time_Data_From_1_D_Dataset(hid_t datasetId, float** DataArray){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::FLOAT_Time_Data_From_1_D_Dataset");
	cout << "DatasetId            " << datasetId << endl;
#endif

	H5Drefresh(datasetId);
	DataSpace = H5Dget_space(datasetId);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
	*DataArray= new float[dims1_out[0]];
	offset1[0]=0;   	MemSpace = H5Screate_simple(1,dims1_out,NULL);
	H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset1,NULL,dims1_out,NULL);
	H5Dread(datasetId, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, *DataArray); 
	H5Sclose(MemSpace); status=H5Sclose(DataSpace);

}



void pvESSI::Interpolate_Stress_Field_At_Nodes(int Time_Index1, int Time_Index2, float Interpolation_Func1,float Interpolation_Func2, float** DataArray){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Interpolate_Stress_Field_At_Nodes");
#endif

	if(Time_Index1==Time_Index2 or Interpolation_Func1==1){
		 float *DataArray1; Extract_Stress_Field_At_Nodes(Time_Index1, &DataArray1);
		 *DataArray = DataArray1;
		 return;
	}
	else if(Interpolation_Func2==1){
		float *DataArray2; Extract_Stress_Field_At_Nodes(Time_Index2, &DataArray2);
		*DataArray = DataArray2;
		return;
	}

	float *DataArray1; Extract_Stress_Field_At_Nodes(Time_Index1, &DataArray1);
	float *DataArray2; Extract_Stress_Field_At_Nodes(Time_Index2, &DataArray2);

	*DataArray = new float[dims2_out[0]];
	for (int i=0;i<dims2_out[0];i++){
		(*DataArray)[i]=Interpolation_Func1*DataArray1[i]+Interpolation_Func2*DataArray2[i];
	}

	delete []  DataArray1; DataArray1 = NULL;
	delete []  DataArray2; DataArray2 = NULL;
}

void pvESSI::Write_Stress_Field_At_Nodes(int TimeIndex1, float **DataArray){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Write_Stress_Field_At_Nodes");
#endif

	H5Fclose(this->id_File);
	this->id_File = H5Fopen(Domain_FileName.c_str(), id_H5F_READ_WRITE, id_H5F_CLOSE_STRONG);

	if(this->id_File>0)
	{

		this->id_Stress_and_Strain = H5Dopen(id_File, "/pvESSI/Field_at_Nodes/Stress_And_Strain", H5P_DEFAULT);
		this->id_Whether_Stress_Strain_Build = H5Dopen(id_File, "/pvESSI/Field_at_Nodes/Whether_Stress_Strain_Build", H5P_DEFAULT);

		offset3[0]   = 0;  					     				 offset3[1]   = TimeIndex1;                             offset3[2] = 0;
	    count3 [0]   = Domain_Number_of_Nodes[this->domain_no];  count3 [2]   = this->Number_of_Strain_Strain_Info;     count3 [1] = 1;
	    dims2_out[0] = Domain_Number_of_Nodes[this->domain_no];	 dims2_out[1] = this->Number_of_Strain_Strain_Info;

	    DataSpace = H5Dget_space(id_Stress_and_Strain);
	    MemSpace = H5Screate_simple(2,dims2_out,NULL);
	    H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset3,NULL,count3,NULL);
	    H5Dwrite(id_Stress_and_Strain, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, (*DataArray)); 
	    H5Sclose(MemSpace); status=H5Sclose(DataSpace);


		int Whether_Stress_Strain_Build = TimeIndex1;
		count1[0]   =1;
		dims1_out[0]=1;
		index_i = (hsize_t)TimeIndex1;

	    HDF5_Write_INT_Array_Data(id_Whether_Stress_Strain_Build,
		                          1,
		                         dims1_out,
		                         &index_i,
		                         NULL,
		                         count1,
		                         NULL,
		                         &Whether_Stress_Strain_Build); // Whether_Stress_Strain_Build_Index 


	    H5Dclose(this->id_Stress_and_Strain);
	    H5Dclose(this->id_Whether_Stress_Strain_Build);

	    cout << "<<<<pvESSI>>>> Storing Stress Field at nodes for future \n" << endl;
	}
	H5Fclose(this->id_File);


	this->id_File = H5Fopen(Domain_FileName.c_str(), id_H5F_READ_ONLY, id_H5F_CLOSE_STRONG);
	this->openDatasetIds();

}


void pvESSI::Extract_Stress_Field_At_Nodes(int Time_Index1, float **Node_Stress_And_Strain_Field){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Extract_Stress_Field_At_Nodes");
#endif

	/*************** Opening the Datasets ***************************/

	float* DataArray;
	hid_t datasetId;

	if(Domain_Read_Status[this->domain_no])
	{

		count1[0]   =1;
		dims1_out[0]=1;
		index_i = (hsize_t)Time_Index1;

		int Whether_Stress_Strain_Build;
		datasetId = H5Dopen(id_File, "pvESSI/Field_at_Nodes/Whether_Stress_Strain_Build", H5P_DEFAULT);
	    HDF5_Read_INT_Array_Data(datasetId,
		                          1,
		                         dims1_out,
		                         &index_i,
		                         NULL,
		                         count1,
		                         NULL,
		                         &Whether_Stress_Strain_Build); // Whether_Stress_Strain_Build
	    H5Dclose(datasetId);

	    if(Whether_Stress_Strain_Build!=-1){
	    	cout << "<<<<pvESSI>>>> Stress-Strains are allready interpolated and stored for this time_step " << this->Time[Time_Index1] <<" s" << endl;
	    	datasetId = H5Dopen(id_File, "pvESSI/Field_at_Nodes/Stress_And_Strain", H5P_DEFAULT);
	    	FLOAT_Time_Data_From_3_D_Dataset(datasetId,Time_Index1,&DataArray);
	    	H5Dclose(datasetId);
	    	(*Node_Stress_And_Strain_Field) = DataArray;
	    	return;
	    }
	}

	Build_Stress_Field_At_Nodes(Time_Index1,&DataArray);
	(*Node_Stress_And_Strain_Field) = DataArray;

	if(Domain_Write_Status[domain_no]){
		Write_Stress_Field_At_Nodes(TimeIndex1, &DataArray);
	}

}

void pvESSI::Build_Stress_Field_At_Nodes(int Time_Index1, float **Node_Stress_And_Strain_Field){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Stress_Field_At_Nodes");
#endif

	count1[0]   =1;
	dims1_out[0]=1;
	index_i = (hsize_t)Time_Index1;

	int GlobalSize = this->Number_of_Strain_Strain_Info * this->Domain_Number_of_Nodes[this->domain_no];
    (*Node_Stress_And_Strain_Field) = new float[GlobalSize];					

	cout << "<<<<pvESSI>>>> Stress-Strains are not interpolated for this time_step " << this->Time[Time_Index1] <<" s" << endl;
	cout << "<<<<pvESSI>>>> I will calculate \n" << endl;

	// initializing the Node_Stress_And_Strain_Field matrix to zero
	for(int i =0; i<GlobalSize;i++)
		(*Node_Stress_And_Strain_Field)[i]=0.0;
	
	///////////////////////////////////////////  Gauss Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	
    float *Gauss_Outputs; FLOAT_Time_Data_From_2_D_Dataset(id_Gauss_Outputs,Time_Index1, &Gauss_Outputs);
	///////////////////////////////////////////// Need to Extend the Dataset ////////////////////////////////////////////////////////////////////////////////////////////
	
	int connectivity_index=0;
	int gauss_output_index=0;
	int nnodes, ngauss, class_tag;  
	int ESSI_ELEMENT_TAG;

	for (int i = 0; i < this->Domain_Number_of_Elements[this->domain_no]; ++i)
	{

		ESSI_ELEMENT_TAG = this->Domain_Element_Map[domain_no][i];

		class_tag  = Domain_Class_Tags[this->domain_no][i];
		nnodes     = NUMBER_OF_NODES(Element_Desc_Array[class_tag]); // Number of element nodes
		ngauss     = NUMBER_OF_GAUSS(Element_Desc_Array[class_tag]); // Number of element gauss nodes

		std::map<int,double**>::iterator it;
		it = Gauss_To_Node_Interpolation_Map.find(class_tag);

		if(ngauss>0){

			if(it != Gauss_To_Node_Interpolation_Map.end()){

				std::vector<int> Nodes_Connectivity_Order = ESSI_to_VTK_Connectivity.find(nnodes)->second;

				////////////////////// Initializing Containers ///////////////////////////////////////
				double **Stress_Strain_At_Nodes        = new double*[nnodes];
			    double **Stress_Strain_At_Gauss_Points = new double*[ngauss];					
				for(int j = 0; j < nnodes; ++j){
			    	Stress_Strain_At_Nodes[j]        = new double[this->Number_of_Strain_Strain_Info];
				}
				for(int j = 0; j < nnodes; ++j){
			    	Stress_Strain_At_Gauss_Points[j] = new double[this->Number_of_Strain_Strain_Info];
				}

				///////////////////// Calculating Stresses at gauss Points /////////////////////////////////

				double  Var_q, Var_p, Var_Plastic_q, Var_Plastic_p;

				for(int j=0; j< nnodes ; j++){

					int index_2=gauss_output_index+j*18;
					for(int k=0; k< 18 ; k++){
						Stress_Strain_At_Gauss_Points[j][k] = Gauss_Outputs[index_2+k];
					}
					// elastic_strain_v   = 1/3*(Gauss_Outputs[index_2+0]+Gauss_Outputs[index_2+1]+Gauss_Outputs[index_2+2]);
					// elastic_strain_dev = sqrt(3/2* (pow(Gauss_Outputs[index_2+0]-Gauss_Outputs[index_2+1],2) +
					// 						   		pow(Gauss_Outputs[index_2+1]-Gauss_Outputs[index_2+2],2) +
					// 						   		pow(Gauss_Outputs[index_2+2]-Gauss_Outputs[index_2+1],2) + 6*
					// 						   		(	pow(Gauss_Outputs[index_2+3],2)+
					// 						   			pow(Gauss_Outputs[index_2+4],2)+
					// 						   			pow(Gauss_Outputs[index_2+5],2))
					// 						  		)
					// 						 );

					Var_Plastic_p   = -1.0*(Gauss_Outputs[index_2+6]+Gauss_Outputs[index_2+7]+Gauss_Outputs[index_2+8]);
					Var_Plastic_q   = sqrt(2.0/9.0* (pow(Gauss_Outputs[index_2+6]-Gauss_Outputs[index_2+1],7) +
								   					 pow(Gauss_Outputs[index_2+7]-Gauss_Outputs[index_2+2],8) +
								   					 pow(Gauss_Outputs[index_2+8]-Gauss_Outputs[index_2+1],6) + 6*
								   					 (	pow(Gauss_Outputs[index_2+9],2)+
								   					 	pow(Gauss_Outputs[index_2+10],2)+
								   					 	pow(Gauss_Outputs[index_2+11],2))
								  					 )
										  );

					Var_p   	   = -1.0/3.0*(Gauss_Outputs[index_2+12]+Gauss_Outputs[index_2+13]+Gauss_Outputs[index_2+14]);
					Var_q 		   = sqrt(1.0/2.0* (pow(Gauss_Outputs[index_2+12]-Gauss_Outputs[index_2+13],2) +
											   		pow(Gauss_Outputs[index_2+13]-Gauss_Outputs[index_2+14],2) +
											   		pow(Gauss_Outputs[index_2+14]-Gauss_Outputs[index_2+13],2) + 6*
											   		(	pow(Gauss_Outputs[index_2+15],2)+
											   			pow(Gauss_Outputs[index_2+16],2)+
											   			pow(Gauss_Outputs[index_2+17],2))
											  		)
										 );

					Stress_Strain_At_Gauss_Points[j][18] = Var_q;
					Stress_Strain_At_Gauss_Points[j][19] = Var_p;
					Stress_Strain_At_Gauss_Points[j][20] = Var_Plastic_q;
					Stress_Strain_At_Gauss_Points[j][21] = Var_Plastic_p;

				}

				vtkMath::MultiplyMatrix	(it->second,Stress_Strain_At_Gauss_Points,nnodes,nnodes,nnodes,this->Number_of_Strain_Strain_Info,Stress_Strain_At_Nodes);



				///////////////////////// Adding the Calculated Stresses at Nodes //////////////////////////
				int node_no=0;
				for(int j=0; j< nnodes ; j++){
					node_no = this->Domain_Connectivity[domain_no][connectivity_index+j];
					// cout << "connectivity_index " << connectivity_index <<" node_no " << node_no  << " class_tag " << class_tag<< endl;
					for(int k=0; k< this->Number_of_Strain_Strain_Info ; k++){
						(*Node_Stress_And_Strain_Field)[node_no*this->Number_of_Strain_Strain_Info+k] = (*Node_Stress_And_Strain_Field)[node_no*this->Number_of_Strain_Strain_Info+k] + Stress_Strain_At_Nodes[j][k] ;
					}
				}

				for(int j = 0; j < nnodes; ++j){
			    	delete [] Stress_Strain_At_Nodes[j];
			    	delete [] Stress_Strain_At_Gauss_Points[j];
				}
				delete [] Stress_Strain_At_Gauss_Points;Stress_Strain_At_Gauss_Points=NULL;
				delete [] Stress_Strain_At_Nodes;Stress_Strain_At_Nodes=NULL;
			}
			else{

				cout << "<<<<pvESSI>>>> Build_Stress_Field_At_Nodes:: Warning!! Gauss_Interpolation to nodes not implemented for element of Class_Tag " << class_tag << endl;
			}

			gauss_output_index = gauss_output_index + ngauss*18;

		}

		connectivity_index = connectivity_index + nnodes;
	}

	// free the allocated arrays 
	delete [] Gauss_Outputs;Gauss_Outputs=NULL;

	// Now take average of contributions from gauss points 
	double epsilon = 1e-6;
	for(int i =0; i<this->Domain_Number_of_Nodes[this->domain_no] ; i++ ){
		for(int j=0; j<this->Number_of_Strain_Strain_Info; j++)
			(*Node_Stress_And_Strain_Field)[i*this->Number_of_Strain_Strain_Info+j] =  (*Node_Stress_And_Strain_Field)[i*this->Number_of_Strain_Strain_Info+j]/((double)this->Domain_Number_of_Gauss_Elements_Shared[this->domain_no][i]+epsilon);
	}


}


/*******************************************************************************
* Interpolating Stress-Strain at Nodes from gauss Points 
********************************************************************************/
void pvESSI::Build_Node_Stress(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh, int Time_Index1, int Time_Index2, float Interpolation_Func1,float Interpolation_Func2){


#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Stress_Field_At_Nodes");
#endif

 	if(Whether_Node_Mesh_Stress_Attributes_Initialized==false)
 	{
		Elastic_Strain = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Elastic_Strain"]);
		Elastic_Strain->Allocate(Total_Number_of_Nodes*6);
		Elastic_Strain->SetNumberOfValues(Total_Number_of_Nodes*6);
		Elastic_Strain->FillComponent(0,0);Elastic_Strain->FillComponent(3,0);
		Elastic_Strain->FillComponent(1,0);Elastic_Strain->FillComponent(4,0);
		Elastic_Strain->FillComponent(2,0);Elastic_Strain->FillComponent(5,0);


		Plastic_Strain = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Plastic_Strain"]);
		Plastic_Strain->Allocate(Total_Number_of_Nodes*6);
		Plastic_Strain->SetNumberOfValues(Total_Number_of_Nodes*6);
		Plastic_Strain->FillComponent(0,0);Plastic_Strain->FillComponent(3,0);
		Plastic_Strain->FillComponent(1,0);Plastic_Strain->FillComponent(4,0);
		Plastic_Strain->FillComponent(2,0);Plastic_Strain->FillComponent(5,0);

		Stress = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Stress"]);
		Stress->Allocate(Total_Number_of_Nodes*6);
		Stress->SetNumberOfValues(Total_Number_of_Nodes*6);
		Stress->FillComponent(0,0);Stress->FillComponent(3,0);
		Stress->FillComponent(1,0);Stress->FillComponent(4,0);
		Stress->FillComponent(2,0);Stress->FillComponent(5,0);

		q = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["q"]);
		q->Allocate(Total_Number_of_Nodes);
		q->SetNumberOfValues(Total_Number_of_Nodes);
		q->FillComponent(0,0);

		p = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["p"]);
		p->Allocate(Total_Number_of_Nodes);
		p->SetNumberOfValues(Total_Number_of_Nodes);
		p->FillComponent(0,0);

		Plastic_Strain_q = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Plastic_Strain_q"]);
		Plastic_Strain_q->Allocate(Total_Number_of_Nodes);
		Plastic_Strain_q->SetNumberOfValues(Total_Number_of_Nodes);
		Plastic_Strain_q->FillComponent(0,0);

		Plastic_Strain_p = vtkSmartPointer<vtkFloatArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Plastic_Strain_p"]);
		Plastic_Strain_p->Allocate(Total_Number_of_Nodes);
		Plastic_Strain_p->SetNumberOfValues(Total_Number_of_Nodes);
		Plastic_Strain_p->FillComponent(0,0);

		this->Whether_Node_Mesh_Stress_Attributes_Initialized=true;
 	}

 	float *Node_Stress_And_Strain_Field; Interpolate_Stress_Field_At_Nodes(Time_Index1, Time_Index2, Interpolation_Func1,Interpolation_Func2, &Node_Stress_And_Strain_Field);

    float Var_q, Var_p, Var_Plastic_q, Var_Plastic_p;
    int  ESSI_NODE_TAG;

	for(int i=0; i< this->Domain_Number_of_Nodes[domain_no]; i++){	

		ESSI_NODE_TAG = this->Domain_Node_Map[this->domain_no][i];

		float El_Strain_Tuple[6] ={
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+0],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+3],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+4],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+1],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+5],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+2]
		};

		float Pl_Strain_Tuple[6] ={
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+6],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+9],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+10],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+7],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+11],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+8]
		};

		float Stress_Tuple[6] ={
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+12],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+15],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+16],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+13],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+17],
			Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+14]
		};

		Elastic_Strain->InsertTypedTuple (this->ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],El_Strain_Tuple);
		Plastic_Strain->InsertTypedTuple (this->ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],Pl_Strain_Tuple);
		Stress->InsertTypedTuple (this->ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],Stress_Tuple);

		q->InsertValue(this->ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+18]);
		p->InsertValue(this->ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+19]);
		Plastic_Strain_q->InsertValue(this->ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+20]);
		Plastic_Strain_p->InsertValue(this->ESSI_To_VTK_Node_Map[ESSI_NODE_TAG],Node_Stress_And_Strain_Field[i*this->Number_of_Strain_Strain_Info+21]);

	}

	Node_Mesh->GetPointData()->AddArray(Elastic_Strain);
	Node_Mesh->GetPointData()->AddArray(Plastic_Strain);
	Node_Mesh->GetPointData()->AddArray(Stress);

	Node_Mesh->GetPointData()->AddArray(q);
	Node_Mesh->GetPointData()->AddArray(p);
	Node_Mesh->GetPointData()->AddArray(Plastic_Strain_q);
	Node_Mesh->GetPointData()->AddArray(Plastic_Strain_p);

	cout << "<<<<pvESSI>>>> Build_Stress_Field_At_Nodes:: Calculation done for the step no  " << Node_Mesh_Current_Time << endl;

	delete [] Node_Stress_And_Strain_Field; Node_Stress_And_Strain_Field = NULL;

	return;
}


void pvESSI::Build_Shared_Info_Per_Mode(){

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::Build_Shared_Info_Per_Mode");
#endif

	for (int i = 0 ; i<Number_of_Nodes; i++){
	   Number_of_Elements_Shared[i]       = 0;
	   Number_of_gauss_Elements_Shared[i] = 0;
	}

}



/**********************************************************************
* Gets the source filename 
**********************************************************************/
std::string pvESSI::GetSourceFile(std::string filename) {

#ifdef DEBUG_MODE
	PRINT_FUNCTION("pvESSI::GetSourceFile");
#endif

  std::string Source_File="";

  char str[filename.size()+1];//as 1 char space for null is also required
  strcpy(str, filename.c_str());

  char * pch;
  // printf ("Splitting string \"%s\" into tokens:\n",str);
  pch = strtok (str,".");
  while (pch != NULL)
  {
    // printf ("%s\n",pch);
    Source_File = Source_File +std::string(pch)+".";
    if(strcmp(pch, "h5") == 0)
		break;
	// pch = strtok (NULL, " ,.-");
	pch = strtok (NULL, ".");
  }
 
 // cout << "Source_File " <<  Source_File << endl;
  return Source_File;
}


void pvESSI::SetPhysicalElementGroupArrayStatus(const char * name, int status){
        if (status){
                this->Physical_Element_Group->EnableArray(name);
        }
        else{
                this->Physical_Element_Group->DisableArray(name);
        }

        this->Modified();
}

const char* pvESSI::GetPhysicalElementGroupArrayName(int index){
	return this->Physical_Element_Group->GetArrayName(index);
}

int pvESSI::GetNumberOfPhysicalElementGroupArrays(){
	return this->Physical_Element_Group->GetNumberOfArrays();
}

int pvESSI::GetPhysicalElementGroupArrayStatus(const char* name){

	return this->Physical_Element_Group->GetArraySetting(name);

}


void pvESSI::SetPhysicalNodeGroupArrayStatus(const char * name, int status){
        if (status){
                this->Physical_Node_Group->EnableArray(name);
        }
        else{
                this->Physical_Node_Group->DisableArray(name);
        }

        this->Modified();
}

const char* pvESSI::GetPhysicalNodeGroupArrayName(int index){
	return this->Physical_Node_Group->GetArrayName(index);
}

int pvESSI::GetNumberOfPhysicalNodeGroupArrays(){
	return this->Physical_Node_Group->GetNumberOfArrays();
}

int pvESSI::GetPhysicalNodeGroupArrayStatus(const char* name){

	return this->Physical_Node_Group->GetArraySetting(name);

}


/*************** Node Visualization Options ****************************/
void pvESSI::SetVisualizationOptionsOnNodeArrayStatus(const char* name, int status){

}
const char* pvESSI::GetVisualizationOptionsOnNodeArrayName(int index){
	return "sa";
}
int pvESSI::GetNumberOfVisualizationOptionsOnNodeArrays(){
	return 0;
}
int pvESSI::GetVisualizationOptionsOnNodeArrayStatus(const char* name){
	return 0;
}




/*************** Node Visualization Options ****************************/
void pvESSI::SetVisualizationOptionsOnElementArrayStatus(const char* name, int status){

}
const char* pvESSI::GetVisualizationOptionsOnElementArrayName(int index){
	return "sa";
}
int pvESSI::GetNumberOfVisualizationOptionsOnElementArrays(){
	return 0;
}
int pvESSI::GetVisualizationOptionsOnElementArrayStatus(const char* name){
	return 0;
}




/*************** Node Visualization Options ****************************/
void pvESSI::SetVisualizationOptionsOnGaussArrayStatus(const char* name, int status){

}
const char* pvESSI::GetVisualizationOptionsOnGaussArrayName(int index){
	return "sa";
}
int pvESSI::GetNumberOfVisualizationOptionsOnGaussArrays(){
	return 0;
}
int pvESSI::GetVisualizationOptionsOnGaussArrayStatus(const char* name){
	return 0;
}




herr_t pvESSI::op_func (hid_t loc_id, const char *name, const H5O_info_t *info,void *operator_data)
{
    // printf ("/");               /* Print root group in object path */

    /*
     * Check if the current object is the root group, and if not print
     * the full path name and type.
     */
    if (name[0] == '.')         /* Root group, do not print '.' */
        // printf ("  (Group)\n");
        printf ("<<<<pvESSI>>>>  Discovering Physical Groups ........   \n");
    else
        switch (info->type) {
            case H5O_TYPE_GROUP:
                // printf ("%s  (Group)\n", name);
                break;
            case H5O_TYPE_DATASET:
                // printf ("%s  (Dataset)\n", name);
            	printf ("%s  Found \n", name);
                Physical_Group_Container.push_back(name);
                break;
            case H5O_TYPE_NAMED_DATATYPE:
                // printf ("%s  (Datatype)\n", name);
                break;
            default:
                printf ("%s  (Unknown)\n", name);
        }

    return 0;
}
