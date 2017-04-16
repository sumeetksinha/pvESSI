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
#include <sstream>
#include "hdf5.h"
#include "vtkProbeFilter.h"
#include <vtkDelaunay3D.h>
#include "vtkSelectionNode.h"
#include "vtkExtractSelection.h"
#include "vtkSelection.h"
#include "vtkGaussianSplatter.h"

// LINK_LIBRARIES(hdf5_cpp hdf5 )

/************************************************************************************************************************************************/
// cmake .. -DParaView_DIR=~/Softwares/Paraview/Paraview-Build/ -DGHOST_BUILD_CDAWEB=OFF
// cmake .. -DParaView_DIR="/home/sumeet/Softwares/ParaView-v4.4.0" -DGHOST_BUILD_CDAWEB=OFF
/************************************************************************************************************************************************/

#define NUMBER_OF_NODES(class_tag_desc)             (class_tag_desc/(int)pow(10.0,(9-(ELE_TAG_DESC_ENCODING/(int)pow(10.0,7))%100%10)))%((int)pow(10.0,((ELE_TAG_DESC_ENCODING/(int)pow(10.0,7))%100%10   - (ELE_TAG_DESC_ENCODING/(int)pow(10.0,8))%100%10  + 1)))
#define NUMBER_OF_GAUSS(class_tag_desc)             (class_tag_desc/(int)pow(10.0,(9-(ELE_TAG_DESC_ENCODING/(int)pow(10.0,5))%100%10)))%((int)pow(10.0,((ELE_TAG_DESC_ENCODING/(int)pow(10.0,5))%100%10   - (ELE_TAG_DESC_ENCODING/(int)pow(10.0,6))%100%10  + 1))) 
#define NUMBER_OF_ELEMENT_OUTPUTS(class_tag_desc)   (class_tag_desc/(int)pow(10.0,(9-(ELE_TAG_DESC_ENCODING/(int)pow(10.0,3))%100%10)))%((int)pow(10.0,((ELE_TAG_DESC_ENCODING/(int)pow(10.0,3))%100%10   - (ELE_TAG_DESC_ENCODING/(int)pow(10.0,4))%100%10  + 1))) 

std::vector<std::string> pvESSI::Physical_Group_Container;
vtkStandardNewMacro(pvESSI);

pvESSI::pvESSI(){ 

	this->FileName = NULL;
	this->eigen_mode_on = false;
	this->Enable_Initialization_Flag=true;
	this->Whether_Physical_Group_Info_build=false;
	this->SetNumberOfInputPorts(0);
	this->SetNumberOfOutputPorts(1);
	this->set_VTK_To_ESSI_Elements_Connectivity();
	this->Build_Meta_Array_Map();
	this->Physical_Element_Group=vtkSmartPointer<vtkDataArraySelection>::New();
	this->Physical_Node_Group=vtkSmartPointer<vtkDataArraySelection>::New();
	Build_Inverse_Matrices();
	Build_Gauss_To_Node_Interpolation_Map();
}
/*****************************************************************************
* This method responds to the request made by vtk Pipeleine. 
* This method is invoked when the time stamp is changed from paraview VCR.
*****************************************************************************/
int pvESSI::RequestData(vtkInformation *vtkNotUsed(request),vtkInformationVector **vtkNotUsed(inputVector),	vtkInformationVector *outputVector){

	vtkInformation *Node_Mesh = outputVector->GetInformationObject(0);
	// outInfo->Print(std::cout);

  	// cout << "this->Node_Mesh_Current_Time " << Node_Mesh->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()) << endl;

	// Get Current Time;
  	this->Node_Mesh_Current_Time = Time_Map.find( Node_Mesh->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))->second;


	piece_no = Node_Mesh->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
	num_of_pieces = Node_Mesh->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

	// cout << "piece_no " << piece_no <<endl;
	// cout << "num_of_pieces " << num_of_pieces << endl;
	// cout << "Node_Mesh_Current_Time " << Node_Mesh_Current_Time <<endl;

	int start =-1;
	int end   =0;
	if(single_file_visualization_mode){
		if(piece_no>0)
			return 1;
	}
	else{
		start = piece_no*ceil(((double)Number_of_Processes_Used)/((double)num_of_pieces));
		end   = start +ceil(((double)Number_of_Processes_Used)/((double)num_of_pieces));
		end = end > (Number_of_Processes_Used-1) ? Number_of_Processes_Used-1:end; 
	}

	// cout << "Start " << start <<endl;
	// cout << "End "   << end << endl;

	for (int i = start; i<end; i++){

		this->domain_no = i;
		Domain_Initializer(domain_no);

		// cout << "domain_no " << domain_no << endl;

		if (!Whether_Node_Mesh_build[domain_no] and (!Show_Gauss_Mesh_Flag or eigen_mode_on or Enable_Displacement_Probing_Flag)){
			UGrid_Node_Mesh[domain_no]         = vtkSmartPointer<vtkUnstructuredGrid>::New();
			UGrid_Current_Node_Mesh[domain_no] = vtkSmartPointer<vtkUnstructuredGrid>::New();
			this->Get_Node_Mesh(UGrid_Node_Mesh[domain_no]);
		}
		else if (!Whether_Gauss_Mesh_build[domain_no] and Show_Gauss_Mesh_Flag){
			UGrid_Gauss_Mesh[domain_no]         = vtkSmartPointer<vtkUnstructuredGrid>::New();
			UGrid_Current_Node_Mesh[domain_no]  = vtkSmartPointer<vtkUnstructuredGrid>::New();
			this->Get_Gauss_Mesh(UGrid_Gauss_Mesh[domain_no]);
		}
				
		if (!Show_Gauss_Mesh_Flag or eigen_mode_on){
			UGrid_Current_Node_Mesh[domain_no]->ShallowCopy(UGrid_Node_Mesh[domain_no]);
		}
		if (Show_Gauss_Mesh_Flag){
			UGrid_Current_Node_Mesh[domain_no]->ShallowCopy(UGrid_Gauss_Mesh[domain_no]);
		}
	  	// int Clength = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	  	// double* Csteps = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

		// ///////////////////////////////  Printing For Debugging ////////////////////////////////
		// cout << "Number of Nodes "   << " " << Number_of_Nodes << endl;
		// cout << "Pseudo_Number_of_Nodes "  << " " << Pseudo_Number_of_Nodes << endl;
		// cout << "Number of Gauss Points " << Number_of_Gauss_Nodes << endl;
		// cout << "Number of Elements "   << " " << Number_of_Elements << endl;
		// cout << "Pseudo_Number_of_Elements"  << " " << Pseudo_Number_of_Elements << endl;
		// ///////////////////////////////////////////////////////////////////////////////////////

		if(eigen_mode_on){
			Build_Eigen_Modes_Node_Attributes(UGrid_Current_Node_Mesh[domain_no], this->Node_Mesh_Current_Time );	
		}

		else{
			if(!Show_Gauss_Mesh_Flag){
				Build_Node_Attributes(UGrid_Current_Node_Mesh[domain_no], this->Node_Mesh_Current_Time );
				if(Enable_Gauss_To_Node_Interpolation_Flag) Build_Stress_Field_At_Nodes(UGrid_Current_Node_Mesh[domain_no], this->Node_Mesh_Current_Time);
				if(Enable_Physical_Node_Group_Selection_Flag or Enable_Physical_Element_Group_Selection_Flag) Build_Physical_Element_Group_Mesh(UGrid_Current_Node_Mesh[domain_no]);
			}
			else if(Show_Gauss_Mesh_Flag){
				Build_Gauss_Attributes(UGrid_Current_Node_Mesh[domain_no], this->Node_Mesh_Current_Time );
			}
		}

	}

	vtkSmartPointer<vtkUnstructuredGrid>  Chunk_Node_mesh =  vtkSmartPointer<vtkUnstructuredGrid>::New();;

	if(single_file_visualization_mode){
		Chunk_Node_mesh->ShallowCopy(UGrid_Current_Node_Mesh[0]);
	}
	else{
		Merge_Mesh(start,end, Chunk_Node_mesh);
	}

	// get the ouptut pointer to paraview 
	vtkUnstructuredGrid *Output_Node_Mesh = vtkUnstructuredGrid::SafeDownCast(Node_Mesh->Get(vtkDataObject::DATA_OBJECT()));

	Output_Node_Mesh->ShallowCopy(Chunk_Node_mesh);


	return 1;
}

// ****************************************************************************
// * This method is called only once for information about time stamp and extent. 
// * This is the method which is called firt after the default constructor is 
// * initialized
// ****************************************************************************

int pvESSI::RequestInformation( vtkInformation *request, vtkInformationVector **vtkNotUsed(inVec), vtkInformationVector* outVec){

	if(this->Enable_Initialization_Flag) this->Initialize(); this->Enable_Initialization_Flag=false;

	vtkInformation* Node_Mesh = outVec->GetInformationObject(0);

	double Time_range[2]={Time[0],Time[Number_of_Time_Steps-1]};

	Node_Mesh->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),Time, this->Number_of_Time_Steps);
	Node_Mesh->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),Time_range,2);
	Node_Mesh->Set(CAN_HANDLE_PIECE_REQUEST(), 1);

	return 1;
}


void pvESSI::Merge_Mesh(int start, int end, vtkSmartPointer<vtkUnstructuredGrid> Mesh){

  vtkSmartPointer<vtkAppendFilter> AppendFilter =  vtkSmartPointer<vtkAppendFilter>::New();

  vtkAppendFilter::GlobalWarningDisplayOff();

  for (int i = start; i<end; i++){
	  AppendFilter->AddInputData(UGrid_Current_Node_Mesh[i]);
	 // vtkIndent indent;
	 // delaunay3D->PrintSelf(std::cout, indent); 
  }

  AppendFilter->Update();
  Mesh->ShallowCopy( AppendFilter->GetOutput());

  return;

}


/*****************************************************************************
* This method creates and process the request from vtk Pipeleine
*****************************************************************************/
int pvESSI::ProcessRequest(vtkInformation  *request_type, vtkInformationVector  **inInfo, vtkInformationVector *outInfo){

	//////////////////////////////////// Printing for Debugging /////////////////////////////////////////
	// cout << "-----------------------------------------------------------------------------------\n";
	// vtkIndent indent;
 	// request_type->PrintSelf(std::cout, indent);
 	////////////////////////////////////////////////////////////////////////////////////////////////////

  return this->Superclass::ProcessRequest(request_type, inInfo, outInfo);
}

/*****************************************************************************
* This method prints about reader plugin i.e itself
*****************************************************************************/
void pvESSI::PrintSelf(ostream& os, vtkIndent indent){

	this->Superclass::PrintSelf(os,indent);
	os << indent << "File Name: " << (this->FileName ? this->FileName : "(none)") << "\n";
	return;
}

/*****************************************************************************
* This method creates a map for the connectivity order from ESSI connectivity
* to vtk elements connectivity. The connectivity order is stored in\
* ESSI_to_VTK_Element map whose key is the number of connectivity nodes. 
*
* !!!I think a better key is needed so that it can be robust.
*****************************************************************************/
void pvESSI::set_VTK_To_ESSI_Elements_Connectivity(){

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

/*****************************************************************************
* This Method builds node attributes for the current time and pushes it to 
* the <vtkUnstructuredGrid> input vtkobject.
*****************************************************************************/
void pvESSI::Build_Node_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh, int Current_Time){

	////////////////////////////////////////////////// For Reference Displacement  ///////////////////////////////////////////////////////////////////////////
	// Finding out the refernce index for displacement 
	int reference_node_mesh_time = Reference_Displacement_Index_Flag;
	if(Reference_Displacement_Index_Flag>=Number_of_Time_Steps) reference_node_mesh_time=Number_of_Time_Steps;
	else if(Reference_Displacement_Index_Flag<=0) reference_node_mesh_time=0;
	float *Reference_Node_Generalized_Displacements;

	////////////////////////////////////////////////////// Reading Node Attributes /////////////////////////////////////////////////////////////////////////////

	int Number_of_DOFs[Number_of_Nodes];
	H5Dread(id_Number_of_DOFs, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Number_of_DOFs); 

	///////////////////////////////////////////  Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	

	DataSpace = H5Dget_space(id_Generalized_Displacements);
	H5Sget_simple_extent_dims(DataSpace, dims2_out, NULL);	
	float Node_Generalized_Displacements[dims2_out[0]];	
	offset2[0]=0; 					  count2[0] = dims2_out[0];		dims1_out[0]=dims2_out[0];
	offset2[1]=Current_Time; 		  count2[1] = 1;				MemSpace = H5Screate_simple(1,dims1_out,NULL);
	H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset2,NULL,count2,NULL);
	H5Dread(id_Generalized_Displacements, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, Node_Generalized_Displacements); 
	H5Sclose(MemSpace);
	H5Sclose(DataSpace); 

	if(Enable_Relative_Displacement_Flag)
	{
		DataSpace = H5Dget_space(id_Generalized_Displacements);
		H5Sget_simple_extent_dims(DataSpace, dims2_out, NULL);	
		offset2[0]=0; 					              count2[0] = dims2_out[0];		dims1_out[0]=dims2_out[0];
		offset2[1]=reference_node_mesh_time; 		  count2[1] = 1;				MemSpace = H5Screate_simple(1,dims1_out,NULL);
		Reference_Node_Generalized_Displacements = new float[dims2_out[0]];
		H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset2,NULL,count2,NULL);
		H5Dread(id_Generalized_Displacements, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, Reference_Node_Generalized_Displacements); 
		H5Sclose(MemSpace);
		H5Sclose(DataSpace); 
	}

	/////////////////////////////////////////////////////////////// DataSets Visulization at Nodes //////////////////////////////////////////////////////////////////////////////////////
 	
 	Generalized_Displacements = vtkSmartPointer<vtkFloatArray>::New(); 
	this->Set_Meta_Array (Meta_Array_Map["Generalized_Displacements"]);
	Generalized_Displacements->Allocate(Number_of_Nodes*3);
	Generalized_Displacements->SetNumberOfValues(Number_of_Nodes*3);

 	if(Enable_uPU_Visualization_Flag)
 	{
	 	Fluid_Displacements = vtkSmartPointer<vtkFloatArray>::New(); 
		this->Set_Meta_Array(Meta_Array_Map["Fluid_Displacements"]);
		Fluid_Displacements->Allocate(Number_of_Nodes*3);
		Fluid_Displacements->SetNumberOfValues(Number_of_Nodes*3);

		// initializing all components to zero 
		Fluid_Displacements->FillComponent (0, 0);
		Fluid_Displacements->FillComponent (1, 0);
		Fluid_Displacements->FillComponent (2, 0);

	 	Pore_Pressure = vtkSmartPointer<vtkFloatArray>::New(); 
		this->Set_Meta_Array (Meta_Array_Map["Pore_Pressure"]);
		Pore_Pressure->Allocate(Number_of_Nodes);
		Pore_Pressure->SetNumberOfValues(Number_of_Nodes);

		// initializing all components to zero 
		Pore_Pressure->FillComponent (0, 0);
	}


	int index_to_generalized_displacement=0;
	
	if(Enable_Relative_Displacement_Flag==false)
	{
		for (int i = 0; i < Number_of_Nodes; i++){

			float disp[3]={
				Node_Generalized_Displacements[index_to_generalized_displacement  ],
				Node_Generalized_Displacements[index_to_generalized_displacement+1],
				Node_Generalized_Displacements[index_to_generalized_displacement+2]
			};


			Generalized_Displacements->InsertTypedTuple(i,disp);

			// ************************************************************** Displacement, Acceleration and Velocity -> Calculation Formulae *********************************************************

			// Displacement(t) = D_t;
			// Acceleration(t) = (1/del.t)^2*{ D_(t+del.t) - 2*D_t + D_(t-del.t)};
			// Acceleration(t) = 1/2/del.t*{ D_(t+del.t) - D_(t-del.t)};

			/***************************************************************** Add Acceleration and Velocity Arrays Here ****************************************/

			// uPU Visualization Mode

		 	if(Enable_uPU_Visualization_Flag and  Number_of_DOFs[i]==7){
		 		
		 		float disp[3]={
					Node_Generalized_Displacements[index_to_generalized_displacement+4],
					Node_Generalized_Displacements[index_to_generalized_displacement+5],
					Node_Generalized_Displacements[index_to_generalized_displacement+6]
				};

				Fluid_Displacements->InsertTypedTuple(i,disp);
				Pore_Pressure->InsertValue(i,Node_Generalized_Displacements[index_to_generalized_displacement+3]);
		 	}


		 	index_to_generalized_displacement = index_to_generalized_displacement + Number_of_DOFs[i];

		}
	}
	else{

		// cout << "reference_node_mesh_time " << reference_node_mesh_time << endl;

		for (int i = 0; i < Number_of_Nodes; i++){

			float disp[3]={
				Node_Generalized_Displacements[index_to_generalized_displacement  ] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement  ],
				Node_Generalized_Displacements[index_to_generalized_displacement+1] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement+1],
				Node_Generalized_Displacements[index_to_generalized_displacement+2] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement+2]
			};

			Generalized_Displacements->InsertTypedTuple(i,disp);


			// uPU Visualization Mode

		 	if(Enable_uPU_Visualization_Flag and Number_of_DOFs[i]==7){		

		 		float disp[3]={
					Node_Generalized_Displacements[index_to_generalized_displacement+4] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement+4],
					Node_Generalized_Displacements[index_to_generalized_displacement+5] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement+5],
					Node_Generalized_Displacements[index_to_generalized_displacement+6] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement+6]
				};

				Fluid_Displacements->InsertTypedTuple(i,disp);
				Pore_Pressure->InsertValue(i,Node_Generalized_Displacements[index_to_generalized_displacement+3] - Reference_Node_Generalized_Displacements[index_to_generalized_displacement+3]);
		 	}

		 	index_to_generalized_displacement = index_to_generalized_displacement + Number_of_DOFs[i];

		}
	}

	Node_Mesh->GetPointData()->AddArray(Generalized_Displacements);


	if(Enable_uPU_Visualization_Flag){	
		Node_Mesh->GetPointData()->AddArray(Fluid_Displacements);
		Node_Mesh->GetPointData()->AddArray(Pore_Pressure);
	}

	/////////////////////////////////////////////////////////////// Support Reactions  //////////////////////////////////////////////////////////////////////////////////////

	if(enable_support_reactions){

		DataSpace = H5Dget_space(id_Support_Reactions);
		H5Sget_simple_extent_dims(DataSpace, dims2_out, NULL);	
		float Reaction_Forces[Number_of_Constrained_Dofs];
		offset2[0]=0; 					  count2[0] = dims2_out[0];		dims1_out[0]=dims2_out[0];
		offset2[1]=Current_Time; 		  count2[1] = 1;				MemSpace = H5Screate_simple(1,dims1_out,NULL);
		H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset2,NULL,count2,NULL);
		H5Dread(id_Support_Reactions, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, Reaction_Forces); 
		H5Sclose(MemSpace);
		H5Sclose(DataSpace); 

		int Constrained_Nodes[Number_of_Constrained_Dofs];
		H5Dread(id_Constrained_Nodes, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Constrained_Nodes); 

		int Constrained_DOFs[Number_of_Constrained_Dofs];
		H5Dread(id_Constrained_DOFs, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Constrained_DOFs); 

		Support_Reactions = vtkSmartPointer<vtkFloatArray>::New(); 
		this->Set_Meta_Array (Meta_Array_Map["Support_Reactions"]);
		Support_Reactions->Allocate(Number_of_Nodes*3);
		Support_Reactions->SetNumberOfValues(Number_of_Nodes*3);

		Support_Reactions->FillComponent (0, 0);
		Support_Reactions->FillComponent (1, 0);
		Support_Reactions->FillComponent (2, 0);

		for(int i = 0; i<Number_of_Constrained_Dofs; i++){

			if(Constrained_DOFs[i]<3)
				Support_Reactions->SetComponent(Constrained_Nodes[i],Constrained_DOFs[i],Reaction_Forces[i]);

		}

		Node_Mesh->GetPointData()->AddArray(Support_Reactions);

	}

 	return;
}

void pvESSI::Build_Eigen_Modes_Node_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh, int Current_Time){


	////////////////////////////////////////////////////// Reading Node Attributes /////////////////////////////////////////////////////////////////////////////

	int Number_of_DOFs[Number_of_Nodes];
	H5Dread(id_Number_of_DOFs, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Number_of_DOFs); 

	///////////////////////////////////////////  Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	

	DataSpace = H5Dget_space(id_modes);
	H5Sget_simple_extent_dims(DataSpace, dims2_out, NULL);	
	float Node_Generalized_Displacements[dims2_out[0]];
	offset2[0]=0; 					  count2[0] = dims2_out[0];		dims1_out[0]=dims2_out[0];
	offset2[1]=Current_Time; 		  count2[1] = 1;				MemSpace = H5Screate_simple(1,dims1_out,NULL);
	H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset2,NULL,count2,NULL);
	H5Dread(id_modes, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, Node_Generalized_Displacements); 
	H5Sclose(MemSpace);
	H5Sclose(DataSpace); 

	/////////////////////////////////////////////////////////////// DataSets Visulization at Nodes //////////////////////////////////////////////////////////////////////////////////////
 	
 	Generalized_Displacements = vtkSmartPointer<vtkFloatArray>::New(); 
	this->Set_Meta_Array (Meta_Array_Map["Generalized_Displacements"]);
	Generalized_Displacements->Allocate(Number_of_Nodes*3);
	Generalized_Displacements->SetNumberOfValues(Number_of_Nodes*3);

	int index_to_generalized_displacement=0;
	

	for (int i = 0; i < Number_of_Nodes; i++){

		float disp[3]={
			Node_Generalized_Displacements[index_to_generalized_displacement  ],
			Node_Generalized_Displacements[index_to_generalized_displacement+1],
			Node_Generalized_Displacements[index_to_generalized_displacement+2]
		};


		Generalized_Displacements->InsertTypedTuple(i,disp);

		index_to_generalized_displacement = index_to_generalized_displacement + Number_of_DOFs[i];

	}

	Node_Mesh->GetPointData()->AddArray(Generalized_Displacements);

 	return;
}


/*****************************************************************************
* This Method builds gauss attributes for the current time and pushes it to 
* the <vtkUnstructuredGrid> input vtkobject.
*****************************************************************************/
void pvESSI::Build_Gauss_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Gauss_Mesh, int Current_Time){

  	if(Enable_Displacement_Probing_Flag){
  		Build_ProbeFilter_Gauss_Mesh(Gauss_Mesh,Current_Time);
  	}

	////////////////////////////////////////////////////// Reading Element  Attributes /////////////////////////////////////////////////////////////////////////////

    int Element_Class_Tags[Number_of_Elements]; 
	H5Dread(id_Class_Tags, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Class_Tags);

	///////////////////////////////////////////  Gauss Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	
	DataSpace = H5Dget_space(id_Gauss_Outputs);
	H5Sget_simple_extent_dims(DataSpace, dims2_out, NULL);
	float *Gauss_Outputs; Gauss_Outputs = new float[dims2_out[0]];
	offset2[0]=0; 					    count2[0] = dims2_out[0];		dims1_out[0]=dims2_out[0];
	offset2[1]=Current_Time;            count2[1] = 1;					MemSpace = H5Screate_simple(1,dims1_out,NULL);
	H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset2,NULL,count2,NULL);
	H5Dread(id_Gauss_Outputs, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, Gauss_Outputs); 
	H5Sclose(MemSpace); status=H5Sclose(DataSpace);

	Elastic_Strain = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Elastic_Strain"]);

	Plastic_Strain = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Plastic_Strain"]);

	Stress = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Stress"]);

	q = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["q"]);

	p = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["p"]);

	Plastic_Strain_q = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Plastic_Strain_q"]);

	Plastic_Strain_p = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Plastic_Strain_p"]);


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

	for (int i = 0; i < this->Number_of_Elements; i++){

    	ngauss = NUMBER_OF_GAUSS(Element_Desc_Array[Element_Class_Tags[i]]);

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

			Elastic_Strain->InsertTypedTuple (index_to_gauss_output,El_Strain_Tuple);
			Plastic_Strain->InsertTypedTuple (index_to_gauss_output,Pl_Strain_Tuple);
			Stress->InsertTypedTuple (index_to_gauss_output,Stress_Tuple);
			Plastic_Strain_p->InsertValue(index_to_gauss_output,Var_Plastic_p);
			Plastic_Strain_q->InsertValue(index_to_gauss_output,Var_Plastic_q);
			p->InsertValue(index_to_gauss_output,Var_p);	
			q->InsertValue(index_to_gauss_output,Var_q);		
			index_to_gauss_output+=1;

		}

	}

	Gauss_Mesh->GetPointData()->AddArray(Stress);
  	Gauss_Mesh->GetPointData()->AddArray(Plastic_Strain);
  	Gauss_Mesh->GetPointData()->AddArray(Elastic_Strain);
	Gauss_Mesh->GetCellData()->AddArray(Plastic_Strain_p);
  	Gauss_Mesh->GetCellData()->AddArray(Plastic_Strain_q);
  	Gauss_Mesh->GetCellData()->AddArray(p);
  	Gauss_Mesh->GetCellData()->AddArray(q);

	return;
}

/*****************************************************************************
* This Method probes node mesh variables at gauss mesh 
* By Default [prob_type =0], probes only generalized displacement
* But the user can choose an option [prob_type = 1] to probe all variables
*****************************************************************************/
void pvESSI::Build_ProbeFilter_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Probe_Input, int Current_Time){

	vtkSmartPointer<vtkUnstructuredGrid> Probe_Source = vtkSmartPointer<vtkUnstructuredGrid>::New();
	Probe_Source->ShallowCopy(this->UGrid_Node_Mesh[domain_no]);

	////////////////////////////////////// For Debugging ////////////////////////////////
	Probe_Source->GetCellData()->RemoveArray("Material_Tag");
	Probe_Source->GetPointData()->RemoveArray("Node_Tag");
	Probe_Source->GetCellData()->RemoveArray("Element_Tag");
	Probe_Source->GetPointData()->RemoveArray("Boundary_Conditions");
	Probe_Source->GetCellData()->RemoveArray("Class_Tag");
	/////////////////////////////////////////////////////////////////////////////////////

	Build_Node_Attributes(Probe_Source, Current_Time);

	/************* Initializing Probe filter ******************************************/

	vtkSmartPointer<vtkProbeFilter> ProbeFilter = vtkSmartPointer<vtkProbeFilter>::New();
	ProbeFilter->SetSourceData(Probe_Source);
	ProbeFilter->SetInputData(Probe_Input);
	ProbeFilter->Update();

	Probe_Input->ShallowCopy( ProbeFilter->GetOutput());

	  // vtkSmartPointer<vtkGaussianSplatter> splatter = 
	  //   vtkSmartPointer<vtkGaussianSplatter>::New();
	  // splatter->SetInputData(Probe_Input);
	  // splatter->SetSampleDimensions(50,50,50);
	  // splatter->SetRadius(0.5);
	  // splatter->ScalarWarpingOff();
	  // splatter->Update();


	// Probe_Input->ShallowCopy(splatter->GetOutput());

  	return;
}


/*****************************************************************************
* This function creates a Node Mesh i.e a mesh from the given node data 
* and element data from hdf5 file. The function stores the mesh in the 
* given input <vtkUnstructuredGrid> input object.
*****************************************************************************/
void pvESSI::Get_Node_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh){

	//////////////////////////////////////////////////// Reading Map Data ///////////////////////////////////////////////////////////////////////////////////////////

	int Node_Map[Number_of_Nodes];
	H5Dread(id_Node_Map, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Node_Map); 

	int Element_Map[Number_of_Elements];
	H5Dread(id_Element_Map, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Map); 

	/////////////////////////////////////////////////// Reading Element Partitioning Information ///////////////////////////////////////////////////////////////////

	int *Element_Partition;
	int *Inverse_Element_Map;

	if(Number_of_Processes_Used>1)
	{

		DataSpace = H5Dget_space(id_Inverse_Element_Map);
		H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
		Inverse_Element_Map = new int[dims1_out[0]];
		offset1[0]=0;   	MemSpace = H5Screate_simple(1,dims1_out,NULL);
		H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset1,NULL,dims1_out,NULL);
		H5Dread(id_Inverse_Element_Map, H5T_NATIVE_INT, MemSpace, DataSpace, H5P_DEFAULT, Inverse_Element_Map); 
		H5Sclose(MemSpace); status=H5Sclose(DataSpace); H5Dclose(id_Inverse_Element_Map);


		// getting the master file containing physical group
		std::string Source_File = GetSourceFile(this->FileName)+"feioutput";
		hid_t id_Source_File = H5Fopen(Source_File.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		id_Element_Partition = H5Dopen(id_Source_File, "/Model/Elements/Partition", H5P_DEFAULT);

		DataSpace = H5Dget_space(id_Element_Partition);
		H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
		Element_Partition = new int[dims1_out[0]]; 
		offset1[0]=0;   	MemSpace = H5Screate_simple(1,dims1_out,NULL);
		H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset1,NULL,dims1_out,NULL);
		H5Dread(id_Element_Partition, H5T_NATIVE_INT, MemSpace, DataSpace, H5P_DEFAULT, Element_Partition); 
		H5Sclose(MemSpace); status=H5Sclose(DataSpace); H5Dclose(id_Element_Partition);

	}

	////////////////////////////////////////////////////// Reading Node Attributes /////////////////////////////////////////////////////////////////////////////

	DataSpace = H5Dget_space(id_Coordinates);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);	
	float Node_Coordinates[dims1_out[0]];
	H5Sclose(DataSpace); 
	H5Dread(id_Coordinates, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Node_Coordinates); 

	int Constrained_Nodes[Number_of_Constrained_Dofs];
	H5Dread(id_Constrained_Nodes, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Constrained_Nodes); 

	int Constrained_DOFs[Number_of_Constrained_Dofs];
	H5Dread(id_Constrained_DOFs, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Constrained_DOFs); 

	//////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	int Element_Class_Tags[Number_of_Elements];
	H5Dread(id_Class_Tags, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Class_Tags);

	int Element_Material_Tags[Number_of_Elements];
	H5Dread(id_Material_Tags, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Material_Tags);

	DataSpace = H5Dget_space(id_Connectivity);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);	
	int Element_Connectivity[dims1_out[0]];
	H5Sclose(DataSpace); 
	H5Dread(id_Connectivity, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Connectivity); 

	/////////////////////////////////////////////////// Building Nodes ///////////////////////////////////////////////////////////////////////////////////////////

	Node_Tag = vtkSmartPointer<vtkIntArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Node_Tag"]);
	Node_Tag->Allocate(Number_of_Nodes);

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(this->Number_of_Nodes);

	Node_Mesh->Allocate(this->Number_of_Elements);

	int index_to_coordinates;

	for (int i = 0; i < Number_of_Nodes; i++){

		Node_Tag->InsertValue(i,Node_Map[i]);

		points->InsertPoint(i,  Node_Coordinates[3*i  ],
								Node_Coordinates[3*i+1],
								Node_Coordinates[3*i+2]
		);
	}


	Node_Mesh->SetPoints(points);
	Node_Mesh->GetPointData()->AddArray(Node_Tag);

	/////////////////////////////////////////////////// Building Boundary Conditions  ///////////////////////////////////////////////////////////////////////////

	Boundary_Conditions = vtkSmartPointer<vtkIntArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Boundary_Conditions"]);
	Boundary_Conditions->Allocate(Number_of_Nodes*3);
	Boundary_Conditions->SetNumberOfValues(Number_of_Nodes*3);
	Boundary_Conditions->FillComponent (0, 0);
	Boundary_Conditions->FillComponent (1, 0);
	Boundary_Conditions->FillComponent (2, 0);

	for(int i = 0; i<Number_of_Constrained_Dofs; i++){

		if(Constrained_DOFs[i]<3)
			Boundary_Conditions->SetComponent(Constrained_Nodes[i],Constrained_DOFs[i],1);

	}

	Node_Mesh->GetPointData()->AddArray(Boundary_Conditions);

	// ///////////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	Material_Tag = vtkSmartPointer<vtkIntArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Material_Tag"]);
	Material_Tag->Allocate(Number_of_Elements);

	Element_Tag = vtkSmartPointer<vtkIntArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Element_Tag"]);
	Element_Tag->Allocate(Number_of_Elements);

	Class_Tag = vtkSmartPointer<vtkIntArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Class_Tag"]);
	Class_Tag->Allocate(Number_of_Elements);


	if(Number_of_Processes_Used>1){
		Partition_Info = vtkSmartPointer<vtkIntArray> ::New();
		this->Set_Meta_Array (Meta_Array_Map["Partition_Info"]);
		Partition_Info->Allocate(Number_of_Elements);		
	}


	int connectivity_index =0;
	int nnodes = 0;
	int Cell_Type;

	for (int i = 0; i < Number_of_Elements; i++){

        nnodes     = NUMBER_OF_NODES(Element_Desc_Array[Element_Class_Tags[i]]);  // Number of element nodes

		Material_Tag->InsertValue(i,Element_Material_Tags[i]);
		Element_Tag->InsertValue(i,Element_Map[i]);
		Class_Tag->InsertValue(i,Element_Class_Tags[i]);

		if(Number_of_Processes_Used>1) {Partition_Info->InsertValue(i,Element_Partition[Element_Map[i]]); }

		vtkIdType Vertices[nnodes];
		Cell_Type = ESSI_to_VTK_Element.find(nnodes)->second;
		std::vector<int> Nodes_Connectivity_Order = ESSI_to_VTK_Connectivity.find(nnodes)->second;

		for(int j=0; j<nnodes ; j++){
			Vertices[j] = Element_Connectivity[connectivity_index+Nodes_Connectivity_Order[j]];
		}

		Node_Mesh->InsertNextCell(Cell_Type, nnodes, Vertices);

		connectivity_index = connectivity_index + nnodes;

	}
	Node_Mesh->GetCellData()->AddArray(Material_Tag);
	Node_Mesh->GetCellData()->AddArray(Element_Tag);
	Node_Mesh->GetCellData()->AddArray(Class_Tag);
	if(Number_of_Processes_Used>1) Node_Mesh->GetCellData()->AddArray(Partition_Info);

	Whether_Node_Mesh_build[domain_no] = true;

	return;
	
}

/*****************************************************************************
* Builds a mesh of a physical group pr groups
* Uses VTKSelection filter to generate the mesh 
*****************************************************************************/
void  pvESSI::Build_Physical_Element_Group_Mesh(vtkSmartPointer<vtkUnstructuredGrid> NodeMesh){

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
				DataSpace = H5Dget_space(id_Physical_Group_Id);
				H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
				int Individual_Physical_group[dims1_out[0]];
				offset1[0]=0;   	MemSpace = H5Screate_simple(1,dims1_out,NULL); 
				H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset1,NULL,dims1_out,NULL);
				H5Dread(id_Physical_Group_Id, H5T_NATIVE_INT, MemSpace, DataSpace, H5P_DEFAULT, Individual_Physical_group); 
				H5Sclose(MemSpace); status=H5Sclose(DataSpace); H5Dclose(id_Physical_Group_Id);

				for(int l=0; l<dims1_out[0];l++){
					Element_Ids->InsertNextValue(Individual_Physical_group[l]);
					number_of_Element_Ids++;
				}

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
				DataSpace = H5Dget_space(id_Physical_Group_Id);
				H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
				int Individual_Physical_group[dims1_out[0]];
				offset1[0]=0;   	MemSpace = H5Screate_simple(1,dims1_out,NULL); 
				H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset1,NULL,dims1_out,NULL);
				H5Dread(id_Physical_Group_Id, H5T_NATIVE_INT, MemSpace, DataSpace, H5P_DEFAULT, Individual_Physical_group); 
				H5Sclose(MemSpace); status=H5Sclose(DataSpace); H5Dclose(id_Physical_Group_Id);

				for(int l=0; l<dims1_out[0];l++){
					Node_Ids->InsertNextValue(Individual_Physical_group[l]);
					number_of_Node_Ids++;
				}

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

/*****************************************************************************
* This Function builds a gauss mesh from the given gauss mesh coordinates  
* Uses Delaunay3D filter to generate the mesh 
*****************************************************************************/
void pvESSI::Get_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Gauss_Mesh){

	/////////////////////////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	int Element_Class_Tags[Number_of_Elements];
	H5Dread(id_Class_Tags, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Class_Tags);

	DataSpace = H5Dget_space(id_Gauss_Point_Coordinates);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);	
	float Element_Gauss_Point_Coordinates[dims1_out[0]];
	H5Dread(id_Gauss_Point_Coordinates, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Gauss_Point_Coordinates); 
	status=H5Sclose(DataSpace); 

	////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(Number_of_Gauss_Points);
	vtkIdType onevertex;

	int gauss_point_no=0;
	int ngauss=0;

	for (int i=0; i < Number_of_Elements; i++){

		ngauss     = (Element_Desc_Array[Element_Class_Tags[i]]%100000)/100;  // Number of gauss nodes

		for(int j=0; j<ngauss ; j++){
			points->InsertPoint(gauss_point_no, 
								Element_Gauss_Point_Coordinates[3*gauss_point_no  ],
								Element_Gauss_Point_Coordinates[3*gauss_point_no+1],
								Element_Gauss_Point_Coordinates[3*gauss_point_no+2]
							   );

			/// Adding 1D point to mesh 
			onevertex = gauss_point_no;
			Gauss_Mesh->InsertNextCell(VTK_VERTEX, 1, &onevertex);

			// updating the gauss point no.
			gauss_point_no++;
		}

	}

	Gauss_Mesh->SetPoints(points);

	// this->Build_Delaunay3D_Gauss_Mesh(Gauss_Mesh);

	Whether_Gauss_Mesh_build[domain_no]=true;

	return;
	
}

/*****************************************************************************
* Generate a tetrahedral mesh from the given input points, and stores the 
* newly generated mesh in the same input vtkUnstructuredGrid object.  
*****************************************************************************/
void pvESSI::Build_Delaunay3D_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> GaussMesh){

  vtkSmartPointer<vtkDelaunay3D> delaunay3D = vtkSmartPointer<vtkDelaunay3D>::New();

  vtkDelaunay3D::GlobalWarningDisplayOff();

  delaunay3D->SetInputData (GaussMesh);
  delaunay3D->Update();

  // vtkIndent indent;
  // delaunay3D->PrintSelf(std::cout, indent);

  GaussMesh->ShallowCopy( delaunay3D->GetOutput());

  return;

}


/*****************************************************************************
* This function generates the metadatay attributes array all at one place  
* I am thinking it to remove it. Was meant to be public accessibe to all but
* currently takes lot of parameters. 
*****************************************************************************/
void pvESSI::Set_Meta_Array( int Meta_Array_Id ){

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

/*****************************************************************************
* Initializes some variables to default  which would be same across all parts
* of the model.
*****************************************************************************/
void pvESSI::Initialize(){

	/************************* To make Debug on *******************/
	// this->DebugOn();

	/*************************************************************/

	/*********Initializing some parameters of visualization ******/
	this->Build_Map_Status = 0;
	this->Number_of_Strain_Strain_Info = 22;
	this->single_file_visualization_mode = false;

    /***************** File_id **********************************/
    this->id_File = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);;  

	/***************** Time Steps *******************************/
    this->id_Number_of_Time_Steps = H5Dopen(id_File, "/Number_of_Time_Steps", H5P_DEFAULT);   
	H5Dread(id_Number_of_Time_Steps, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Time_Steps);
	H5Dclose(id_Number_of_Time_Steps);

	H5Eset_auto (NULL, NULL, NULL);  // To stop HDF5 from printing error message
	this->id_Eigen_Mode_Analysis                 = H5Gopen(id_File, "Eigen_Mode_Analysis", H5P_DEFAULT);
	if(this->id_Eigen_Mode_Analysis>0){
		this->id_number_of_modes				   = H5Dopen(id_File,"Eigen_Mode_Analysis/number_of_modes",H5P_DEFAULT);
		H5Dread(id_number_of_modes, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Time_Steps);
		H5Dclose(id_number_of_modes);

		cout << "<<<<pvESSI>>>> Eigen_Mode_Analysis is On!!! \n" << endl;
		eigen_mode_on = true;
	}

    this->Time = new double[Number_of_Time_Steps];
    float temp_Time [Number_of_Time_Steps];
	this->id_time = H5Dopen(id_File, "/time", H5P_DEFAULT); 
	H5Dread(id_time, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,temp_Time);
	H5Dclose(id_time);

    /***************** Model Info *******************************/
	this->id_Number_of_Iterations      = H5Dopen(id_File, "/Number_of_Iterations", H5P_DEFAULT);
	H5Dread(id_Number_of_Iterations, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Sub_Steps );
	H5Dclose(id_Number_of_Iterations);

    this->id_Number_of_Processes_Used = H5Dopen(id_File, "/Number_of_Processes_Used", H5P_DEFAULT);
	H5Dread(id_Number_of_Processes_Used, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Processes_Used);
	H5Dclose(id_Number_of_Processes_Used);

    this->id_Process_Number           = H5Dopen(id_File, "/Process_Number", H5P_DEFAULT);
	H5Dread(id_Process_Number, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Process_Number);
	H5Dclose(id_Process_Number);

	/****************** Element Class Description *************************/
	
	// Finiding Number of Elements in ESSI  
    this->id_Element_Class_Desc = H5Dopen(id_File, "Model/Elements/Element_Class_Desc", H5P_DEFAULT);
    DataSpace = H5Dget_space(id_Element_Class_Desc);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
	H5Sclose(DataSpace);

	Element_Desc_Array = new int[dims1_out[0]];
	H5Dread(id_Element_Class_Desc, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Desc_Array);
	H5Dclose(id_Element_Class_Desc);

	ELE_TAG_DESC_ENCODING = Element_Desc_Array[0];  // Class_DESC_ENCODING

	if(ELE_TAG_DESC_ENCODING==-1)
		ELE_TAG_DESC_ENCODING = 235789180;


	H5Fclose(this->id_File);		// close the file

	// Initializing Time vector
	for(int p=0; p<Number_of_Time_Steps;p++)
		this->Time[p]=temp_Time[p];
		// this->Time[p] = p;

	Build_Time_Map(); 

	if(Number_of_Processes_Used==1 || Process_Number>0){
		this->single_file_visualization_mode = true;
		this->Number_of_Processes_Used =1;
	}

	UGrid_Current_Node_Mesh   = new vtkSmartPointer<vtkUnstructuredGrid>[Number_of_Processes_Used];		   // Holds the current combination of node mesh
	// UGrid_Current_Gauss_Mesh  = new vtkSmartPointer<vtkUnstructuredGrid>[Number_of_Processes_Used];		   // Hold the current combination of gauss mesh
	UGrid_Node_Mesh = new vtkSmartPointer<vtkUnstructuredGrid>[Number_of_Processes_Used];
	UGrid_Gauss_Mesh = new vtkSmartPointer<vtkUnstructuredGrid>[Number_of_Processes_Used];
	Whether_Node_Mesh_build = new bool[Number_of_Processes_Used];
	Whether_Gauss_Mesh_build = new bool[Number_of_Processes_Used];

	for(int i=0; i<Number_of_Processes_Used; i++ ){
		Whether_Node_Mesh_build[i] = false;
		Whether_Gauss_Mesh_build[i]= false;
	}


}

void pvESSI::Domain_Initializer(int Domain_Number){

	std::string filename;

	if(Domain_Number>=0){
		std::string Source_File = GetSourceFile(this->FileName);
		std::stringstream ss;
		int digits = Number_of_Processes_Used > 0 ? (int) log10 ((double) Number_of_Processes_Used) + 1 : 1;
		ss << setfill('0') << setw(digits) << Domain_Number+1;
		filename = Source_File + ss.str()+".feioutput";
	}
	else{
		filename = this->FileName; 
		domain_no =0;               // setting the index to be zero for Single_Domain_Visualization_Mode
	}

	if(!Whether_Physical_Group_Info_build)
	{
		//********************* Building Avialable Physical Groups ***************************************************************/
		// getting the master file containing physical group
		std::string Source_File = GetSourceFile(this->FileName)+"feioutput";
		hid_t id_Source_File = H5Fopen(Source_File.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	    this->id_Physical_Element_Groups = H5Gopen(id_Source_File, "/Model/Physical_Groups/Physical_Element_Groups", H5P_DEFAULT);
	    this->id_Physical_Node_Groups    = H5Gopen(id_Source_File, "/Model/Physical_Groups/Physical_Node_Groups", H5P_DEFAULT);

	    status = H5Ovisit (id_Physical_Element_Groups, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, NULL);
		for(int i=0; i<Physical_Group_Container.size(); ++i){ const char *name = Physical_Group_Container[i].c_str(); this->Physical_Element_Group->AddArray(name);}
		if(!Enable_Physical_Element_Group_Selection_Flag){this->Physical_Element_Group->DisableAllArrays();} Physical_Group_Container.clear();

	    status = H5Ovisit (id_Physical_Node_Groups, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, NULL);
	    for(int i=0; i<Physical_Group_Container.size(); ++i){ const char *name = Physical_Group_Container[i].c_str(); this->Physical_Node_Group->AddArray(name);}
	    if(!Enable_Physical_Node_Group_Selection_Flag){this->Physical_Node_Group->DisableAllArrays();} Physical_Group_Container.clear();

	    Whether_Physical_Group_Info_build = true;;

	    H5Fclose(id_Source_File);
		H5Gclose(id_Physical_Element_Groups); 
		H5Gclose(id_Physical_Node_Groups); 
	}

	/****************************************************************************************
	* We want to open the file and keep it open until we have read all the data what we want
	* So lets go wohoo!!!
	****************************************************************************************/

	  /***************** File_id **********************************/
	  cout << filename.c_str() <<endl;
	  this->id_File = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	  /***************** Read Model Information *******************/
      this->id_Number_of_Elements       = H5Dopen(id_File, "/Number_of_Elements", H5P_DEFAULT); 
      H5Dread(id_Number_of_Elements, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Elements );
      H5Dclose(id_Number_of_Elements);
    
      this->id_Number_of_Nodes          = H5Dopen(id_File, "/Number_of_Nodes", H5P_DEFAULT);
      H5Dread(id_Number_of_Nodes, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Nodes );
      H5Dclose(id_Number_of_Nodes);
    
      this->id_Number_of_Gauss_Points   = H5Dopen(id_File, "/Number_of_Gauss_Points", H5P_DEFAULT);
      H5Dread(id_Number_of_Gauss_Points, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Number_of_Gauss_Points );
      H5Dclose(id_Number_of_Gauss_Points);

	  /****************** Check if pvESSI build *******************/
	  H5Eset_auto (NULL, NULL, NULL);  // To stop HDF% from printing error message
	  this->id_pvESSI = H5Gopen(id_File, "/pvESSI", H5P_DEFAULT);  

	  if(id_pvESSI>0 and Enable_Building_of_Maps_Flag) { // Every thing is ready to go

        cout << "<<<<pvESSI>>>>  Maps are build HURRAY!!! \n" << endl;

		// Finiding Number of connectivity nodes
		this->id_Connectivity = H5Dopen(id_File, "Model/Elements/Connectivity", H5P_DEFAULT);
		DataSpace = H5Dget_space(id_Connectivity);
		H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
		this->Number_of_Connectivity_Nodes = dims1_out[0];
		H5Sclose(DataSpace);

		// Finding Number of Constrained nodes
		this->id_Constrained_Nodes = H5Dopen(id_File, "Model/Nodes/Constrained_Nodes", H5P_DEFAULT);
		DataSpace = H5Dget_space(id_Constrained_Nodes);
		H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
		this->Number_of_Constrained_Dofs = dims1_out[0];
		H5Sclose(DataSpace);
		H5Dclose(id_Constrained_Nodes);

	  }
	  else{ // We need to build the pvESSI folder

		cout << "<<<<pvESSI>>>>  Maps Not Build, So lets go and build it first \n " << endl;
		H5Fclose(this->id_File);
		this->id_File = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		Build_Maps();
		H5Fclose(this->id_File);
		this->id_File = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		Enable_Building_of_Maps_Flag = true;
		cout << "<<<<pvESSI>>>>  Maps are now built \n " << endl;

	  }

	  if(eigen_mode_on){
	  	this->id_frequencies					   = H5Dopen(id_File,"Eigen_Mode_Analysis/frequencies",H5P_DEFAULT); 	
		this->id_modes 							   = H5Dopen(id_File,"Eigen_Mode_Analysis/modes",H5P_DEFAULT);
		this->id_periods						   = H5Dopen(id_File,"Eigen_Mode_Analysis/periods",H5P_DEFAULT);
		this->id_values  						   = H5Dopen(id_File,"Eigen_Mode_Analysis/values",H5P_DEFAULT);
	  }

	  /**************** Open Datasets in pvESSI Folder *************************************/
		  this->id_Element_Map                     = H5Dopen(id_File, "pvESSI/Element_Map", H5P_DEFAULT);
		  this->id_Connectivity                    = H5Dopen(id_File, "pvESSI/Connectivity", H5P_DEFAULT);
		  this->id_Node_Map                        = H5Dopen(id_File, "pvESSI/Node_Map", H5P_DEFAULT);
		  this->id_Number_of_DOFs                  = H5Dopen(id_File, "pvESSI/Number_of_DOFs", H5P_DEFAULT);
		  this->id_Inverse_Node_Map                = H5Dopen(id_File, "pvESSI/Inverse_Node_Map", H5P_DEFAULT);
		  this->id_Inverse_Element_Map             = H5Dopen(id_File, "pvESSI/Inverse_Element_Map", H5P_DEFAULT);
		  this->id_Number_of_Elements_Shared       = H5Dopen(id_File, "pvESSI/Number_of_Elements_Shared", H5P_DEFAULT);
		  this->id_Number_of_Gauss_Elements_Shared = H5Dopen(id_File, "pvESSI/Number_of_Gauss_Elements_Shared", H5P_DEFAULT);
		  this->id_Constrained_Nodes               = H5Dopen(id_File, "pvESSI/Constrained_Nodes", H5P_DEFAULT);
		  this->id_Class_Tags                      = H5Dopen(id_File, "pvESSI/Class_Tags", H5P_DEFAULT);
		  this->id_Material_Tags                   = H5Dopen(id_File, "pvESSI/Material_Tags", H5P_DEFAULT);
	      this->id_Physical_Element_Groups         = H5Gopen(id_File, "pvESSI/Physical_Groups/Physical_Element_Groups", H5P_DEFAULT);
	      this->id_Physical_Node_Groups            = H5Gopen(id_File, "pvESSI/Physical_Groups/Physical_Node_Groups", H5P_DEFAULT);


		  /*************** Field at Nodes ***************************/
		  this->id_Field_at_Nodes_group = H5Gopen(id_File, "pvESSI/Field_at_Nodes", H5P_DEFAULT);
		  this->id_Stress_and_Strain = H5Dopen(id_File, "pvESSI/Field_at_Nodes/Stress_And_Strain", H5P_DEFAULT);
		  this->id_Whether_Stress_Strain_Build = H5Dopen(id_File, "pvESSI/Field_at_Nodes/Whether_Stress_Strain_Build", H5P_DEFAULT);
		  // this->id_Energy = H5Dopen(id_File, "/Field_at_Nodes/Energy", H5P_DEFAULT);                                         // Not implemented
		  // this->id_Whether_Energy_Build = H5Dopen(id_File, "/Field_at_Nodes/Whether_Energy_Build", H5P_DEFAULT);             // Not implemented
		  

	  /**************** Element Info ******************************/
	  this->id_Gauss_Point_Coordinates = H5Dopen(id_File, "Model/Elements/Gauss_Point_Coordinates", H5P_DEFAULT);
	  // this->id_Element_Outputs = H5Dopen(id_File, "Model/Elements/Element_Outputs", H5P_DEFAULT);
	  this->id_Gauss_Outputs = H5Dopen(id_File, "Model/Elements/Gauss_Outputs", H5P_DEFAULT);
	  // this->id_Substep_Outputs = H5Dopen(id_File, "Model/Elements/Substep_Outputs", H5P_DEFAULT);

	  /**************** Node Info ******************************/
	  this->id_Constrained_DOFs = H5Dopen(id_File, "Model/Nodes/Constrained_DOFs", H5P_DEFAULT);
	  this->id_Coordinates = H5Dopen(id_File, "Model/Nodes/Coordinates", H5P_DEFAULT);
	  this->id_Generalized_Displacements = H5Dopen(id_File, "Model/Nodes/Generalized_Displacements", H5P_DEFAULT);
	  this->id_Support_Reactions = H5Dopen(id_File,"/Model/Nodes/Support_Reactions", H5P_DEFAULT);
	  if(id_Support_Reactions>0) enable_support_reactions=true; else enable_support_reactions=false;

}


// pvESSI:: ~pvESSI(){
// 	this->Close_File();
// }

void pvESSI::Close_File(){

	  /***************** File_id **********************************/
	  H5Fclose(this->id_File);

   //   /**************** Element Info ******************************/
	  // H5Gclose(this->id_Elements_group); 
	  // H5Dclose(this->id_Class_Tags);
	  // H5Dclose(this->id_Connectivity);
	  // H5Dclose(this->id_Element_types);
	  // H5Dclose(this->id_Gauss_Point_Coordinates);
	  // H5Dclose(this->id_Index_to_Connectivity);
	  // H5Dclose(this->id_Index_to_Gauss_Point_Coordinates);
	  // H5Dclose(this->id_Index_to_Outputs);
	  // H5Dclose(this->id_Material_tags);
	  // H5Dclose(this->id_Number_of_Gauss_Points);
	  // H5Dclose(this->id_Element_Number_of_Nodes);
	  // H5Dclose(this->id_Number_of_Output_Fields);
	  // H5Dclose(this->id_Outputs);
	  // // H5Dclose(this->id_Substep_Output);

	  // /**************** Node Info ******************************/
	  // H5Gclose(this->id_Nodes_group); 
	  // H5Dclose(this->id_Constrained_DOFs);
	  // H5Dclose(this->id_Constarined_Nodes);
	  // H5Dclose(this->id_Coordinates);
	  // H5Dclose(this->id_Generalized_Displacements);
	  // H5Dclose(this->id_Index_to_Coordinates);
	  // H5Dclose(this->id_Index_to_Generalized_Displacements);
	  // H5Dclose(this->id_Number_of_DOFs);

	  // /**************** Maps ***********************************/
	  // H5Gclose(this->id_Maps_group); 
	  // H5Dclose(this->id_Element_Map);
	  // H5Dclose(this->id_Node_Map);
	  // H5Dclose(this->id_Inverse_Node_Map);
	  // H5Dclose(this->id_Inverse_Element_Map);
	  // H5Dclose(this->id_Number_of_Elements_Shared);
	  // H5Dclose(this->id_Number_of_Gauss_Elements_Shared);

	  // /*************** Field at Nodes ***************************/
	  // H5Gclose(this->id_Field_at_Nodes_group);
	  // H5Dclose(this->id_Stress_and_Strain);
	  // H5Dclose(this->id_Whether_Stress_Strain_Build);
	  // // H5Dclose(this->id_Energy);                           // Not implemented
	  // // H5Dclose(this->id_Whether_Energy_Build);             // Not implemented
	  
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


void pvESSI::Build_Maps(){


	// Finiding Pseudo number of nodes 
    this->id_Number_of_DOFs = H5Dopen(id_File, "Model/Nodes/Number_of_DOFs", H5P_DEFAULT);
    DataSpace = H5Dget_space(id_Number_of_DOFs);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
	this->Pseudo_Number_of_Nodes = dims1_out[0];
	H5Sclose(DataSpace);

	// Finiding Pseudo number of elements 
	this->id_Class_Tags = H5Dopen(id_File, "Model/Elements/Class_Tags", H5P_DEFAULT);
	DataSpace = H5Dget_space(id_Class_Tags);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
	this->Pseudo_Number_of_Elements = dims1_out[0];
	H5Sclose(DataSpace);

	// Finiding Number of connectivity nodes
	this->id_Connectivity = H5Dopen(id_File, "Model/Elements/Connectivity", H5P_DEFAULT);
	DataSpace = H5Dget_space(id_Connectivity);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
	this->Number_of_Connectivity_Nodes = dims1_out[0];
	H5Sclose(DataSpace);

	// Finding Number of Constrained nodes
	this->id_Constrained_Nodes = H5Dopen(id_File, "Model/Nodes/Constrained_Nodes", H5P_DEFAULT);
	DataSpace = H5Dget_space(id_Constrained_Nodes);
	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
	this->Number_of_Constrained_Dofs = dims1_out[0];
	H5Sclose(DataSpace);

    this->id_Material_Tags           = H5Dopen(id_File, "Model/Elements/Material_Tags", H5P_DEFAULT);

	/*********************************************** First Need to build Datasets and Folders *********************************************************************/

	id_pvESSI = H5Gcreate(id_File, "/pvESSI", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
	H5Gclose(id_pvESSI); 
 
	id_Field_at_Nodes_group = H5Gcreate(id_File, "/pvESSI/Field_at_Nodes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
	H5Gclose(id_Field_at_Nodes_group); 

	dims1_out[0]= Number_of_Nodes;

	DataSpace = H5Screate_simple(1, dims1_out, NULL);
	id_Node_Map = H5Dcreate(id_File,"pvESSI/Node_Map",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(DataSpace);
	H5Dclose(id_Node_Map);
	id_Node_Map = H5Dopen(id_File,"pvESSI/Node_Map", H5P_DEFAULT);

	DataSpace = H5Screate_simple(1, dims1_out, NULL);
	hid_t pvESSI_id_Number_of_DOFs = H5Dcreate(id_File,"pvESSI/Number_of_DOFs",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(DataSpace);
	H5Dclose(pvESSI_id_Number_of_DOFs);
	pvESSI_id_Number_of_DOFs = H5Dopen(id_File,"pvESSI/Number_of_DOFs", H5P_DEFAULT);

	DataSpace = H5Screate_simple(1, dims1_out, NULL);
	id_Number_of_Elements_Shared = H5Dcreate(id_File,"pvESSI/Number_of_Elements_Shared",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(DataSpace);
	H5Dclose(id_Number_of_Elements_Shared);
	id_Number_of_Elements_Shared = H5Dopen(id_File,"pvESSI/Number_of_Elements_Shared", H5P_DEFAULT);

	DataSpace = H5Screate_simple(1, dims1_out, NULL);
	id_Number_of_Gauss_Elements_Shared = H5Dcreate(id_File,"pvESSI/Number_of_Gauss_Elements_Shared",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(DataSpace);
	H5Dclose(id_Number_of_Gauss_Elements_Shared);
	id_Number_of_Gauss_Elements_Shared = H5Dopen(id_File,"pvESSI/Number_of_Gauss_Elements_Shared", H5P_DEFAULT);

	dims1_out[0]= Number_of_Elements; 
	DataSpace = H5Screate_simple(1, dims1_out, NULL);
	id_Element_Map = H5Dcreate(id_File,"pvESSI/Element_Map",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(DataSpace);	
	H5Dclose(id_Element_Map);
	id_Element_Map = H5Dopen(id_File,"pvESSI/Element_Map", H5P_DEFAULT);

	DataSpace = H5Screate_simple(1, dims1_out, NULL);
	hid_t pvESSI_id_Class_Tags = H5Dcreate(id_File,"pvESSI/Class_Tags",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(DataSpace);
	H5Dclose(pvESSI_id_Class_Tags);
	pvESSI_id_Class_Tags = H5Dopen(id_File,"pvESSI/Class_Tags", H5P_DEFAULT);

	DataSpace = H5Screate_simple(1, dims1_out, NULL);
	hid_t pvESSI_id_Material_Tags = H5Dcreate(id_File,"pvESSI/Material_Tags",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(DataSpace);
	H5Dclose(pvESSI_id_Material_Tags);
	pvESSI_id_Material_Tags = H5Dopen(id_File,"pvESSI/Material_Tags", H5P_DEFAULT);

	dims1_out[0]= Number_of_Connectivity_Nodes; 
	DataSpace = H5Screate_simple(1, dims1_out, NULL);
	hid_t pvESSI_id_Connectivity = H5Dcreate(id_File,"pvESSI/Connectivity",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(DataSpace);
	H5Dclose(pvESSI_id_Connectivity);
	pvESSI_id_Connectivity = H5Dopen(id_File,"pvESSI/Connectivity", H5P_DEFAULT);

	dims1_out[0]= Pseudo_Number_of_Nodes; 
	DataSpace = H5Screate_simple(1, dims1_out, NULL);
	id_Inverse_Node_Map = H5Dcreate(id_File,"pvESSI/Inverse_Node_Map",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(DataSpace);
	H5Dclose(id_Inverse_Node_Map);
	id_Inverse_Node_Map = H5Dopen(id_File,"pvESSI/Inverse_Node_Map", H5P_DEFAULT);

	dims1_out[0]= Pseudo_Number_of_Elements; 
	DataSpace = H5Screate_simple(1, dims1_out, NULL);
	id_Inverse_Element_Map = H5Dcreate(id_File,"pvESSI/Inverse_Element_Map",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(DataSpace);
	H5Dclose(id_Inverse_Element_Map);
	id_Inverse_Element_Map = H5Dopen(id_File,"pvESSI/Inverse_Element_Map", H5P_DEFAULT);

	dims1_out[0]= Number_of_Constrained_Dofs; 
	DataSpace = H5Screate_simple(1, dims1_out, NULL);
	hid_t pvESSI_id_Constrained_Nodes = H5Dcreate(id_File,"pvESSI/Constrained_Nodes",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(DataSpace);
	H5Dclose(pvESSI_id_Constrained_Nodes);
	pvESSI_id_Constrained_Nodes = H5Dopen(id_File,"pvESSI/Constrained_Nodes", H5P_DEFAULT);

	/************************************************* Write default Values of -1 for all time steps *****************************************/
	dims1_out[0]= this->Number_of_Time_Steps; 
	DataSpace = H5Screate_simple(1, dims1_out, NULL);
	id_Whether_Stress_Strain_Build = H5Dcreate(id_File,"pvESSI/Field_at_Nodes/Whether_Stress_Strain_Build",H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(DataSpace);
	H5Dclose(id_Whether_Stress_Strain_Build);
	id_Whether_Stress_Strain_Build = H5Dopen(id_File,"pvESSI/Field_at_Nodes/Whether_Stress_Strain_Build", H5P_DEFAULT);

	int Whether_Stress_Strain_Build[Number_of_Time_Steps];
	int index = -1;
	while(++index<Number_of_Time_Steps){
		Whether_Stress_Strain_Build[index]=-1;
	}

	count1[0]   =Number_of_Time_Steps;
	dims1_out[0]=Number_of_Time_Steps;
	index_i     =0;

    HDF5_Write_INT_Array_Data(id_Whether_Stress_Strain_Build, 1, dims1_out, &index_i, NULL, count1, NULL, Whether_Stress_Strain_Build); // Whether_Stress_Strain_Build

	// *********************************************************** Creating Strain Dataset ****************************************************
	dims3[0]=this->Number_of_Nodes;       dims3[1]= this->Number_of_Time_Steps;   dims3[2]= Number_of_Strain_Strain_Info; 
	maxdims3[0]=this->Number_of_Nodes; maxdims3[1]= this->Number_of_Time_Steps;   maxdims3[2]= Number_of_Strain_Strain_Info;
	DataSpace = H5Screate_simple(3, dims3, maxdims3);

    hid_t prop = H5Pcreate (H5P_DATASET_CREATE);
    hsize_t chunk_dims[3] = {this->Number_of_Nodes,1,Number_of_Strain_Strain_Info};
    H5Pset_chunk (prop, 3, chunk_dims);
	id_Stress_and_Strain = H5Dcreate(id_File,"pvESSI/Field_at_Nodes/Stress_And_Strain",H5T_NATIVE_DOUBLE,DataSpace,H5P_DEFAULT,prop, H5P_DEFAULT); 
	status = H5Sclose(DataSpace);
	H5Dclose(id_Stress_and_Strain);
	id_Stress_and_Strain = H5Dopen(id_File,"pvESSI/Field_at_Nodes/Stress_And_Strain", H5P_DEFAULT);


	//************** Building Node Map datset *******************************//
	int *Node_Map = new int[Number_of_Nodes];
	int *pvESSI_Number_of_DOFs = new int[Number_of_Nodes];
	int *Inverse_Node_Map = new int[Pseudo_Number_of_Nodes];
	int *pvESSI_Constrained_Nodes = new int[Number_of_Constrained_Dofs];

    int *Number_of_DOFs; Number_of_DOFs  = new int[Pseudo_Number_of_Nodes];
	H5Dread(id_Number_of_DOFs, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Number_of_DOFs);

    int *Constrained_Nodes; Constrained_Nodes = new int[Number_of_Constrained_Dofs];
	H5Dread(id_Constrained_Nodes, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Constrained_Nodes);

	index   = -1;
	node_no = 0;
	while(++index < this->Pseudo_Number_of_Nodes){

	     if(Number_of_DOFs[index]!=-1){

	     	Node_Map[node_no] = index;
	     	pvESSI_Number_of_DOFs[node_no] = Number_of_DOFs[index];
	     	Inverse_Node_Map[index] = node_no++; 

	     }
	     else{

	     	Inverse_Node_Map[index]=-1;
	     }
	}

	for(int i =0; i<Number_of_Constrained_Dofs; i++){
		pvESSI_Constrained_Nodes[i] = Inverse_Node_Map[Constrained_Nodes[i]];
	}

	// write data to the datasets
	H5Dwrite(id_Node_Map, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Node_Map);
	H5Dwrite(id_Inverse_Node_Map, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Inverse_Node_Map);
	H5Dwrite(pvESSI_id_Number_of_DOFs, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,pvESSI_Number_of_DOFs);
	H5Dwrite(pvESSI_id_Constrained_Nodes, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,pvESSI_Constrained_Nodes);

	// close the datasets 
	// H5Dclose(id_Number_of_DOFs);
	// H5Dclose(id_Constrained_Nodes);


	//************** Building Element Map datset *******************************//
	int *Element_Map = new int[Number_of_Elements];
	int *pvESSI_Connectivity = new int[Number_of_Connectivity_Nodes];
	int *pvESSI_Class_Tags = new int[Number_of_Elements];
	int *pvESSI_Material_Tags = new int[Pseudo_Number_of_Elements];
	int *Inverse_Element_Map = new int[Pseudo_Number_of_Elements];
	int *Number_of_Elements_Shared = new int[Number_of_Nodes];
	int *Number_of_gauss_Elements_Shared =  new int[Number_of_Nodes];

	for (int i = 0 ; i<Number_of_Nodes; i++){
		Number_of_Elements_Shared[i]       = 0;
		Number_of_gauss_Elements_Shared[i] = 0;
	}

	int *Element_Class_Tags = new int[Pseudo_Number_of_Elements];
	H5Dread(id_Class_Tags, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Class_Tags);

	int *Element_Connectivity = new int[Number_of_Connectivity_Nodes];
	H5Dread(id_Connectivity, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Connectivity); 

	this->id_Index_to_Connectivity = H5Dopen(id_File, "Model/Elements/Index_to_Connectivity", H5P_DEFAULT);
	int *Element_Index_to_Connectivity = new int[Pseudo_Number_of_Elements];
	H5Dread(id_Index_to_Connectivity, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Index_to_Connectivity); 

	int *Material_Tags = new int[Pseudo_Number_of_Elements];
	H5Dread(id_Material_Tags, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Material_Tags); 

	std::map<int,double**>::iterator it;

	int nnodes     = 0;
	int class_tag  = 0;
	int class_desc = 0;
	int ngauss     = 0;
	int index_to_connectivity = 0 ;

	int index_2;
	index   = -1;
	element_no = 0;

	while(++index < this->Pseudo_Number_of_Elements){

		class_tag  = Element_Class_Tags[index];
		class_desc = Element_Desc_Array[class_tag];
        nnodes     = NUMBER_OF_NODES(class_desc); // Number of element nodes
        ngauss     = NUMBER_OF_GAUSS(class_desc); // Number of gauss nodes

        // cout <<  nnodes << " " << ngauss << endl;

	     if(class_tag!=-1){

	     	Element_Map[element_no] = index;
	     	pvESSI_Class_Tags[element_no] = class_tag;
	     	pvESSI_Material_Tags[element_no] = Material_Tags[index];
	     	Inverse_Element_Map[index] = element_no++; 


	     	// cout << " Class_Tag[index] " << Element_Class_Tags[index] <<endl; 
	     	// cout << " Number of elemnts " << Int_Variable_2 << endl;
	     	// cout << " Index to connectivity " << Element_Index_to_Connectivity[index] << endl;

	     	if(ngauss>0){

	     		it = Gauss_To_Node_Interpolation_Map.find(class_tag);

				   if(it != Gauss_To_Node_Interpolation_Map.end()){

				   		index_2= -1;
						while(++index_2 <nnodes){

							node_no = Inverse_Node_Map[Element_Connectivity[index_2+index_to_connectivity]];
							Number_of_gauss_Elements_Shared[node_no] = Number_of_gauss_Elements_Shared[node_no]+1;
							Number_of_Elements_Shared[node_no]      = Number_of_Elements_Shared[node_no]+1;
							pvESSI_Connectivity[index_2+index_to_connectivity] = node_no;
						}
				   }
				   else{
				   	   	
				   	   	index_2= -1;
						while(++index_2 <nnodes){
							node_no = Inverse_Node_Map[Element_Connectivity[index_2+index_to_connectivity]];
							Number_of_Elements_Shared[node_no]      = Number_of_Elements_Shared[node_no]+1;
							pvESSI_Connectivity[index_2+index_to_connectivity] = node_no;
						}

				   	   cout << "<<<<[pvESSI]>>>> Build_Maps:: Warning!! Gauss to Node Interpolation not implemented for element of Class_Tag  " << Element_Class_Tags[index] << endl;

				   }

	     	}
		     else{

			   	   	index_2= -1;
					while(++index_2 <nnodes){
						node_no = Inverse_Node_Map[Element_Connectivity[index_2+index_to_connectivity]];
						Number_of_Elements_Shared[node_no]      = Number_of_Elements_Shared[node_no]+1;
						pvESSI_Connectivity[index_2+index_to_connectivity] = node_no;
					}

		     }

		    index_to_connectivity = index_to_connectivity + nnodes;
	     }
	     else{

	     	Inverse_Element_Map[index]=-1;
	     }
	}

	H5Dwrite(id_Element_Map, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Map);
	H5Dwrite(id_Inverse_Element_Map, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Inverse_Element_Map);
	H5Dwrite(id_Number_of_Elements_Shared, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Number_of_Elements_Shared);
	H5Dwrite(id_Number_of_Gauss_Elements_Shared, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Number_of_gauss_Elements_Shared);
	H5Dwrite(pvESSI_id_Connectivity, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,pvESSI_Connectivity);
	H5Dwrite(pvESSI_id_Class_Tags, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,pvESSI_Class_Tags);
	H5Dwrite(pvESSI_id_Material_Tags, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,pvESSI_Material_Tags);

    // close dataset 
	H5Dclose(id_Connectivity);
	H5Dclose(id_Class_Tags);
	H5Dclose(id_Material_Tags);

	/**************************************************** Building Element Physical Group *****************************************************/

	// getting the master file containing physical group
	std::string Source_File = GetSourceFile(this->FileName)+"feioutput";
	hid_t id_Source_File = H5Fopen(Source_File.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    this->id_Physical_Element_Groups = H5Gopen(id_Source_File, "/Model/Physical_Groups/Physical_Element_Groups", H5P_DEFAULT);
    this->id_Physical_Node_Groups    = H5Gopen(id_Source_File, "/Model/Physical_Groups/Physical_Node_Groups", H5P_DEFAULT);

	hid_t pvESSI_id_Physical_Groups = H5Gcreate(id_File, "/pvESSI/Physical_Groups", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 	                                        
	hid_t pvESSI_id_Physical_Element_Groups = H5Gcreate(id_File, "/pvESSI/Physical_Groups/Physical_Element_Groups", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);             
	hid_t pvESSI_id_Physical_Node_Groups = H5Gcreate(id_File, "/pvESSI/Physical_Groups/Physical_Node_Groups", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);   

    hid_t id_Physical_Group_Id;
    int Length_of_individual_physical_group=0;

	//**********/ Physical Element Groups **********************//
    Physical_Group_Container.clear();
    status = H5Ovisit (id_Physical_Element_Groups, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, NULL);

    for(int i=0; i<Physical_Group_Container.size(); ++i){

    	const char *name = Physical_Group_Container[i].c_str();

    	id_Physical_Group_Id =  H5Dopen(id_Physical_Element_Groups,name, H5P_DEFAULT);
    	DataSpace = H5Dget_space(id_Physical_Group_Id);
    	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
    	int *Individual_Physical_group = new int[dims1_out[0]];
    	int *pvESSI_Individual_Physical_group = new int[dims1_out[0]]; 
    	offset1[0]=0;   	MemSpace = H5Screate_simple(1,dims1_out,NULL);
		H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset1,NULL,dims1_out,NULL);
		H5Dread(id_Physical_Group_Id, H5T_NATIVE_INT, MemSpace, DataSpace, H5P_DEFAULT, Individual_Physical_group); 
		H5Sclose(MemSpace); status=H5Sclose(DataSpace); H5Dclose(id_Physical_Group_Id);

		Length_of_individual_physical_group =0;
		for(int j=0; j<(int)dims1_out[0];j++){
			element_no = Individual_Physical_group[j];
			if(Element_Class_Tags[element_no]!=-1){
				// cout << Individual_Physical_group[j] << endl;;
				// cout << Inverse_Element_Map[element_no] << endl;
				pvESSI_Individual_Physical_group[Length_of_individual_physical_group++]=element_no;
			}
		}

		dims1_out[0] = Length_of_individual_physical_group;
		DataSpace = H5Screate_simple(1, dims1_out, NULL);
		id_Physical_Group_Id = H5Dcreate(pvESSI_id_Physical_Element_Groups,name,H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(DataSpace);

		if(Length_of_individual_physical_group>0)
		{
			H5Dwrite(id_Physical_Group_Id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,pvESSI_Individual_Physical_group);
		}

		H5Dclose(id_Physical_Group_Id);

		delete [] pvESSI_Individual_Physical_group;
		delete [] Individual_Physical_group;
    }

	//**********/ Physical Node Groups **********************//
    Physical_Group_Container.clear();
    status = H5Ovisit (id_Physical_Node_Groups, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, NULL);

    for(int i=0; i<Physical_Group_Container.size(); ++i){

    	const char *name = Physical_Group_Container[i].c_str();

    	id_Physical_Group_Id =  H5Dopen(id_Physical_Node_Groups,name, H5P_DEFAULT);
    	DataSpace = H5Dget_space(id_Physical_Group_Id);
    	H5Sget_simple_extent_dims(DataSpace, dims1_out, NULL);
    	int *Individual_Physical_group = new int[dims1_out[0]];
    	int *pvESSI_Individual_Physical_group = new int[dims1_out[0]];
    	offset1[0]=0;   	MemSpace = H5Screate_simple(1,dims1_out,NULL);
		H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset1,NULL,dims1_out,NULL);
		H5Dread(id_Physical_Group_Id, H5T_NATIVE_INT, MemSpace, DataSpace, H5P_DEFAULT, Individual_Physical_group); 
		H5Sclose(MemSpace); status=H5Sclose(DataSpace); H5Dclose(id_Physical_Group_Id);

		Length_of_individual_physical_group =0;
		for(int j=0; j<dims1_out[0];j++){
			node_no = Individual_Physical_group[j];
			if(Number_of_DOFs[node_no]!=-1){
				// cout << Individual_Physical_group[j] << endl;;
				// cout << Inverse_Node_Map[node_no] << endl;
				pvESSI_Individual_Physical_group[Length_of_individual_physical_group++]=node_no;
			}
		}

		dims1_out[0] = (hsize_t) Length_of_individual_physical_group;
		DataSpace = H5Screate_simple(1, dims1_out, NULL);
		id_Physical_Group_Id = H5Dcreate(pvESSI_id_Physical_Node_Groups,name,H5T_NATIVE_INT,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
		H5Sclose(DataSpace);


		if(Length_of_individual_physical_group>0)
		{
			H5Dwrite(id_Physical_Group_Id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,pvESSI_Individual_Physical_group);
		}

		H5Dclose(id_Physical_Group_Id);

		delete [] pvESSI_Individual_Physical_group;
		delete [] Individual_Physical_group;
    }

	H5Gclose(pvESSI_id_Physical_Groups); 
	H5Gclose(pvESSI_id_Physical_Element_Groups); 
	H5Gclose(pvESSI_id_Physical_Node_Groups); 
	H5Gclose(id_Physical_Element_Groups); 
	H5Gclose(id_Physical_Node_Groups); 
}

/*****************************************************************************
* Builds time map i.e. a way to get the index number from a given non-interger
* time request from paraview desk. Although it works fine, its not the great 
* way to achieve it. 
* Needs to be changed. I think, paraview has fixed this issue in their latest 
* version, but I am not sure, need to check it. 
*****************************************************************************/
void pvESSI::Build_Time_Map(){

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
* Builds Meta Array map, so that it is easie to add attributes to vtkObject 
* and is neat and clean.
*****************************************************************************/
void pvESSI::Build_Meta_Array_Map(){

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

	free(Temp_Twenty_Seven_Brick);
	free(Temp_Twenty_Node_Brick);
	free(Temp_Eight_Node_Brick);


	return;
  }

void pvESSI::Build_Brick_Coordinates(){

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


void pvESSI::HDF5_Read_FLOAT_Array_Data(hid_t id_DataSet, int rank, hsize_t *data_dims, hsize_t *offset, hsize_t *stride, hsize_t *count, hsize_t *block, float* data){

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
	      H5T_NATIVE_FLOAT,     // Format of data in memory
	      MemSpace,           // Description of data in memory
	      DataSpace,          // Description of data in storage (including selection)
	      H5P_DEFAULT,         // Form of reading
	      data                // The actual data
	  );

  H5Sclose(DataSpace);
  H5Sclose(MemSpace);
}

void pvESSI::HDF5_Write_FLOAT_Array_Data(hid_t id_DataSet, int rank, hsize_t *data_dims, hsize_t *offset, hsize_t *stride, hsize_t *count, hsize_t *block, float* data){

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
	      id_DataSet,        // Dataset to write to
	      H5T_NATIVE_FLOAT,     // Format of data in memory
	      MemSpace,           // Description of data in memory
	      DataSpace,          // Description of data in storage (including selection)
	      H5P_DEFAULT,         // Form of reading
	      data                // The actual data
	  );

  H5Sclose(DataSpace);
  H5Sclose(MemSpace);
}


void pvESSI::HDF5_Read_DOUBLE_Array_Data(hid_t id_DataSet, int rank, hsize_t *data_dims, hsize_t *offset, hsize_t *stride, hsize_t *count, hsize_t *block, double* data){

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
	      id_DataSet,        // Dataset to write to
	      H5T_NATIVE_DOUBLE,     // Format of data in memory
	      MemSpace,           // Description of data in memory
	      DataSpace,          // Description of data in storage (including selection)
	      H5P_DEFAULT,         // Form of reading
	      data                // The actual data
	  );

  H5Sclose(DataSpace);
  H5Sclose(MemSpace);
}


void pvESSI::HDF5_Write_DOUBLE_Array_Data(hid_t id_DataSet, int rank, hsize_t *data_dims, hsize_t *offset, hsize_t *stride, hsize_t *count, hsize_t *block, double* data){

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
	      id_DataSet,        // Dataset to write to
	      H5T_NATIVE_DOUBLE,     // Format of data in memory
	      MemSpace,           // Description of data in memory
	      DataSpace,          // Description of data in storage (including selection)
	      H5P_DEFAULT,         // Form of reading
	      data                // The actual data
	  );

  H5Sclose(DataSpace);
  H5Sclose(MemSpace);
}

/*******************************************************************************
* Interpolating Stress-Strain at Nodes from gauss Points 
********************************************************************************/
void pvESSI::Build_Stress_Field_At_Nodes(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh, int Node_Mesh_Current_Time){

	// First Check whether Stresses has been calculated or not.
	// If Not Calculated : Calculate it and store it.
	// if Calculated: Visualize it

	Elastic_Strain = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Elastic_Strain"]);

	Plastic_Strain = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Plastic_Strain"]);

	Stress = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Stress"]);

	q = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["q"]);

	p = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["p"]);

	Plastic_Strain_q = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Plastic_Strain_q"]);

	Plastic_Strain_p = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Plastic_Strain_p"]);

	count1[0]   =1;
	dims1_out[0]=1;
	index_i = (hsize_t)Node_Mesh_Current_Time;

	int Whether_Stress_Strain_Build;
    HDF5_Read_INT_Array_Data(id_Whether_Stress_Strain_Build,
	                          1,
	                         dims1_out,
	                         &index_i,
	                         NULL,
	                         count1,
	                         NULL,
	                         &Whether_Stress_Strain_Build); // Whether_Stress_Strain_Build

	double Node_Stress_And_Strain_Field[this->Number_of_Nodes][Number_of_Strain_Strain_Info];

	if(Whether_Stress_Strain_Build==-1){

		cout << "<<<<pvESSI>>>> Stress-Strains are not interpolated for this time_step " << Node_Mesh_Current_Time <<endl;
		cout << "<<<<pvESSI>>>> I will calculate and store it for future \n" << endl;

		//////////////////////////////////////////////// Extending the dataset /////////////////////////////////////////////////////////////////////////////
		
		// DataSpace = H5Dget_space (id_Stress_and_Strain);
	 //    H5Sget_simple_extent_dims (DataSpace, dims3, NULL);
		// H5Sclose(DataSpace);

		// initializing the Node_Stress_And_Strain_Field matrix to zero
		for(int i =0; i<this->Number_of_Nodes;i++)
			for(int j =0; j<Number_of_Strain_Strain_Info; j++)
				Node_Stress_And_Strain_Field[i][j]=0.0;
		
	    // dims3[0] = this->Number_of_Nodes;
	    // dims3[1] = dims3[1]<Number_of_Time_Steps?(dims3[1]+1):dims3[1];
	    // dims3[2] = this->Number_of_Strain_Strain_Info;
	    // status = H5Dset_extent (id_Stress_and_Strain, dims3);

		offset3[0] = 0;  					     offset3[1] = Node_Mesh_Current_Time;    offset3[2] = 0;
	    count3 [0] = this->Number_of_Nodes;		 count3 [1] = 1;		                 count3 [2] = this->Number_of_Strain_Strain_Info;
	    dims2_out[0] =this->Number_of_Nodes;	 dims2_out[1] = this->Number_of_Strain_Strain_Info;

	    DataSpace = H5Dget_space(id_Stress_and_Strain);
	    MemSpace = H5Screate_simple(2,dims2_out,NULL);
	    H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset3,NULL,count3,NULL);
	    H5Dwrite(id_Stress_and_Strain, H5T_NATIVE_DOUBLE, MemSpace, DataSpace, H5P_DEFAULT, Node_Stress_And_Strain_Field); 
	    H5Sclose(MemSpace); status=H5Sclose(DataSpace);
		
		////////////////////////////////////////////////// Reading Map Data ///////////////////////////////////////////////////////////////////////////////////////////

		int Number_of_Gauss_Elements_Shared[Number_of_Nodes];
		H5Dread(id_Number_of_Gauss_Elements_Shared, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Number_of_Gauss_Elements_Shared);

		////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	    int Element_Class_Tags[Number_of_Elements]; 
		H5Dread(id_Class_Tags, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Class_Tags);

		int Element_Connectivity[Number_of_Connectivity_Nodes];
		H5Dread(id_Connectivity, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Connectivity); 

		///////////////////////////////////////////  Gauss Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	
		DataSpace = H5Dget_space(id_Gauss_Outputs);
		H5Sget_simple_extent_dims(DataSpace, dims2_out, NULL);
 		float *Gauss_Outputs; Gauss_Outputs =  new float[dims2_out[0]];
		offset2[0]=0; 					    count2[0] = dims2_out[0];		dims1_out[0]=dims2_out[0];
		offset2[1]=Node_Mesh_Current_Time;  count2[1] = 1;					MemSpace = H5Screate_simple(1,dims1_out,NULL);
		H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset2,NULL,count2,NULL);
		H5Dread(id_Gauss_Outputs, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, Gauss_Outputs); 
		H5Sclose(MemSpace); status=H5Sclose(DataSpace);

		// ///////////////////////////////////////////// Need to Extend the Dataset ////////////////////////////////////////////////////////////////////////////////////////////
		
		// for(int i =0 ; i<this->Number_of_Nodes;i++){
		// 	for(int j=0;j<18;j++)
		// 		Node_Stress_And_Strain_Field[i][j]=0;
		// 		// cout <<  Node_Stress_And_Strain_Field[i][j] << " " ;
		// 	cout << endl;
		// }

		int connectivity_index=0;
		int gauss_output_index=0;
		int nnodes, ngauss, class_tag;  

		for (int i = 0; i < this->Number_of_Elements; ++i)
		{

			class_tag  = Element_Class_Tags[i];
			nnodes     = NUMBER_OF_NODES(Element_Desc_Array[class_tag]); // Number of element nodes
			ngauss     = NUMBER_OF_GAUSS(Element_Desc_Array[class_tag]); // Number of element gauss nodes

			// cout << class_tag << "  " << nnodes << " " << ngauss << " " << endl;

			std::map<int,double**>::iterator it;
			it = Gauss_To_Node_Interpolation_Map.find(class_tag);

			if(ngauss>0){

				if(it != Gauss_To_Node_Interpolation_Map.end()){

					std::vector<int> Nodes_Connectivity_Order = ESSI_to_VTK_Connectivity.find(nnodes)->second;

					////////////////////// Initializing Containers ///////////////////////////////////////
					double **Stress_Strain_At_Nodes = new double*[nnodes];
				    double **Stress_Strain_At_Gauss_Points = new double*[ngauss];					
					for(int j = 0; j < nnodes; ++j){
				    	Stress_Strain_At_Nodes[j] = new double[this->Number_of_Strain_Strain_Info];
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

					/*******************************************************************************************/

					// vtkIdType Vertices[No_of_Element_Nodes];
					// Cell_Type = ESSI_to_VTK_Element.find(No_of_Element_Nodes)->second;
					// std::vector<int> Nodes_Connectivity_Order = ESSI_to_VTK_Connectivity.find(No_of_Element_Nodes)->second;

					// for(int j=0; j<No_of_Element_Nodes ; j++){
					// 	Vertices[j] = Inverse_Node_Map[Element_Connectivity[connectivity_index+Nodes_Connectivity_Order[j]]];
					// }

					// *****************************************************************************************

					// for(int p =0; p<number_of_element_nodes ; p++ ){
					// 	for(int q=0; q<18; q++)
					// 		cout << Stress_Strain_At_Gauss_Points[p][q] << " ";
					// 	cout << endl;
					// }
					// cout << endl;

					vtkMath::MultiplyMatrix	(it->second,Stress_Strain_At_Gauss_Points,nnodes,nnodes,nnodes,this->Number_of_Strain_Strain_Info,Stress_Strain_At_Nodes);

					///////////////////////// Adding the Calculated Stresses at Nodes //////////////////////////
					int node_no=0;
					for(int j=0; j< nnodes ; j++){
						node_no = Element_Connectivity[connectivity_index+j];
						// cout << "connectivity_index " << connectivity_index <<" node_no " << node_no  << " class_tag " << class_tag<< endl;
						for(int k=0; k< this->Number_of_Strain_Strain_Info ; k++){
							Node_Stress_And_Strain_Field[node_no][k] = Node_Stress_And_Strain_Field[node_no][k] + Stress_Strain_At_Nodes[j][k] ;
						}
					}

				}
				else{

					cout << "<<<<pvESSI>>>> Build_Stress_Field_At_Nodes:: Warning!! Gauss_Interpolation to nodes not implemented for element of Class_Tag " << class_tag << endl;
				}

				gauss_output_index = gauss_output_index + ngauss*18;

			}
			connectivity_index = connectivity_index + nnodes;
		}

		// free(Gauss_Outputs);	

		// Now take average of contributions from gauss points 
		
		double epsilon = 1e-6;
		for(int i =0; i<this->Number_of_Nodes ; i++ ){
			for(int j=0; j<this->Number_of_Strain_Strain_Info; j++)
				Node_Stress_And_Strain_Field[i][j] =  Node_Stress_And_Strain_Field[i][j]/((double)Number_of_Gauss_Elements_Shared[i]+epsilon);
		}

		offset3[0] = 0;  					     offset3[1] = Node_Mesh_Current_Time;    offset3[2] = 0;
	    count3 [0] = this->Number_of_Nodes;		 count3 [1] = 1;		    count3 [2] = this->Number_of_Strain_Strain_Info;
	    dims2_out[0] =this->Number_of_Nodes;	 dims2_out[1] = this->Number_of_Strain_Strain_Info;

	    DataSpace = H5Dget_space(id_Stress_and_Strain);
	    MemSpace = H5Screate_simple(2,dims2_out,NULL);
	    H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset3,NULL,count3,NULL);
	    H5Dwrite(id_Stress_and_Strain, H5T_NATIVE_DOUBLE, MemSpace, DataSpace, H5P_DEFAULT, Node_Stress_And_Strain_Field); 
	    H5Sclose(MemSpace); status=H5Sclose(DataSpace);


		Whether_Stress_Strain_Build = Node_Mesh_Current_Time;
		count1[0]   =1;
		dims1_out[0]=1;
		index_i = (hsize_t)Node_Mesh_Current_Time;

	    HDF5_Write_INT_Array_Data(id_Whether_Stress_Strain_Build,
		                          1,
		                         dims1_out,
		                         &index_i,
		                         NULL,
		                         count1,
		                         NULL,
		                         &Whether_Stress_Strain_Build); // Whether_Stress_Strain_Build_Index 

	}

	///////////////////// Reading the stress-strain at nodes ///////////////////////////

	offset3[0] = 0;  					     offset3[1] =Whether_Stress_Strain_Build;    offset3[2] = 0;
    count3 [0] = this->Number_of_Nodes;		 count3 [1] = 1;		    	             count3 [2] = this->Number_of_Strain_Strain_Info;
    dims2_out[0] =this->Number_of_Nodes;	 dims2_out[1] = this->Number_of_Strain_Strain_Info;

    DataSpace = H5Dget_space(id_Stress_and_Strain);
    MemSpace = H5Screate_simple(2,dims2_out,NULL);
    H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset3,NULL,count3,NULL);
    H5Dread(id_Stress_and_Strain, H5T_NATIVE_DOUBLE, MemSpace, DataSpace, H5P_DEFAULT, Node_Stress_And_Strain_Field); 
    H5Sclose(MemSpace); status=H5Sclose(DataSpace);

    float Var_q, Var_p, Var_Plastic_q, Var_Plastic_p;


	for(int i=0; i< this->Number_of_Nodes; i++){	

		float El_Strain_Tuple[6] ={
			Node_Stress_And_Strain_Field[i][0],
			Node_Stress_And_Strain_Field[i][3],
			Node_Stress_And_Strain_Field[i][4],
			Node_Stress_And_Strain_Field[i][1],
			Node_Stress_And_Strain_Field[i][5],
			Node_Stress_And_Strain_Field[i][2]
		};

		float Pl_Strain_Tuple[6] ={
			Node_Stress_And_Strain_Field[i][6],
			Node_Stress_And_Strain_Field[i][9],
			Node_Stress_And_Strain_Field[i][10],
			Node_Stress_And_Strain_Field[i][7],
			Node_Stress_And_Strain_Field[i][11],
			Node_Stress_And_Strain_Field[i][8]
		};

		float Stress_Tuple[6] ={
			Node_Stress_And_Strain_Field[i][12],
			Node_Stress_And_Strain_Field[i][15],
			Node_Stress_And_Strain_Field[i][16],
			Node_Stress_And_Strain_Field[i][13],
			Node_Stress_And_Strain_Field[i][17],
			Node_Stress_And_Strain_Field[i][14]
		};

		Elastic_Strain->InsertTypedTuple (i,El_Strain_Tuple);
		Plastic_Strain->InsertTypedTuple (i,Pl_Strain_Tuple);
		Stress->InsertTypedTuple (i,Stress_Tuple);

		q->InsertValue(i,Node_Stress_And_Strain_Field[i][18]);
		p->InsertValue(i,Node_Stress_And_Strain_Field[i][19]);
		Plastic_Strain_q->InsertValue(i,Node_Stress_And_Strain_Field[i][20]);
		Plastic_Strain_p->InsertValue(i,Node_Stress_And_Strain_Field[i][21]);
	}

	Node_Mesh->GetPointData()->AddArray(Elastic_Strain);
	Node_Mesh->GetPointData()->AddArray(Plastic_Strain);
	Node_Mesh->GetPointData()->AddArray(Stress);

	Node_Mesh->GetPointData()->AddArray(q);
	Node_Mesh->GetPointData()->AddArray(p);
	Node_Mesh->GetPointData()->AddArray(Plastic_Strain_q);
	Node_Mesh->GetPointData()->AddArray(Plastic_Strain_p);
	cout << "<<<<pvESSI>>>> Build_Stress_Field_At_Nodes:: Calculation done for the step no  " << Node_Mesh_Current_Time << endl;

	return;
}

/**********************************************************************
* Gets the source filename 
**********************************************************************/
std::string pvESSI::GetSourceFile(std::string filename) {

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
