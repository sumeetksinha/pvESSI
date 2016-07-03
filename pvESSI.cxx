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

using std::ostringstream;


// LINK_LIBRARIES(hdf5_cpp hdf5 )

/************************************************************************************************************************************************/
// cmake .. -DParaView_DIR=~/Softwares/Paraview/Paraview-Build/ -DGHOST_BUILD_CDAWEB=OFF
// cmake .. -DParaView_DIR="/home/sumeet/Softwares/ParaView-v4.4.0" -DGHOST_BUILD_CDAWEB=OFF
/************************************************************************************************************************************************/

using namespace::H5;

vtkStandardNewMacro(pvESSI);
 
pvESSI::pvESSI(){ 

	this->FileName = NULL;
	this->SetNumberOfInputPorts(0);
	this->SetNumberOfOutputPorts(1);
	UGrid_Gauss_Mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
	UGrid_Node_Mesh  = vtkSmartPointer<vtkUnstructuredGrid>::New();
	UGrid_All_Mesh   = vtkSmartPointer<vtkUnstructuredGrid>::New();
	this->set_VTK_To_ESSI_Elements_Connectivity();
	this->initialize();
}

int pvESSI::RequestData(vtkInformation *vtkNotUsed(request),vtkInformationVector **vtkNotUsed(inputVector),	vtkInformationVector *outputVector){
 
	// get the info object
	vtkInformation *Node_Mesh = outputVector->GetInformationObject(0);
	// vtkInformation *Gauss_Mesh = outputVector->GetInformationObject(1);
	// outInfo->Print(std::cout);

	if (!Whether_Node_Mesh_Build){this->Get_Node_Mesh();} 
	
	// if (!Whether_Gauss_Mesh_Build ){ this->Get_Gauss_Mesh();} 

	// int extent[6] = {0,-1,0,-1,0,-1};
	// outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);

  	Node_Mesh_Current_Time = Time_Map.find( Node_Mesh->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))->second;
  	// Gauss_Mesh_Current_Time = Time_Map.find( Node_Mesh->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))->second;

  	// int Clength = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  	// double* Csteps = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

 	// cout<< "Node_Mesh_Current_Time " << Node_Mesh->Get(vtkStreamingDemandDrivenPipeline::DATA_TIME_STEP())<< endl;


 	// cout<< "Gauss_Mesh_Current_Time " << Gauss_Mesh_Current_Time<< endl;
	// cout<< "Clength " << "  " << vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP() << endl;
	// for (int i =0 ; i< Clength ; i++){
	// 	cout << Csteps[i] << "  " << endl;
	// }

    // cout << Node_Mesh_Current_Time/Time[0] << endl;

	// if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
	//          std::cout << "Update Time" << endl;
	//     if (!outInfo->Has(vtkDataObject::DATA_TIME_STEP()))
	    		// std::cout << "DATA_TIME_STEP" << endl;

	// return outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())
	// int piece, numPieces;
	// piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
	// numPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

	/////////////////////////////// Read in Paralll //////////////////////////////
	// outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
	///////////////////////////////////////////////////////////////////////////////

	// get the ouptut pointer to paraview 
	// if(Node_Mesh_Current_Time!=Node_Mesh_Previous_Time){
		Build_Node_Attributes();
	// 	Node_Mesh_Previous_Time = Node_Mesh_Current_Time;
	// }
	// else if(Gauss_Mesh_Current_Time!=Gauss_Mesh_Previous_Time){
		// Build_Gauss_Attributes();
	// 	Gauss_Mesh_Previous_Time = Gauss_Mesh_Current_Time;
	// }

	vtkUnstructuredGrid *UGrid_Current_Node_Mesh = vtkUnstructuredGrid::SafeDownCast(Node_Mesh->Get(vtkDataObject::DATA_OBJECT()));
	// vtkUnstructuredGrid *UGrid_Current_Gauss_Mesh = vtkUnstructuredGrid::SafeDownCast(Gauss_Mesh->Get(vtkDataObject::DATA_OBJECT()));
	UGrid_Current_Node_Mesh->ShallowCopy(UGrid_Node_Mesh);
	// UGrid_Current_Gauss_Mesh->ShallowCopy(UGrid_Gauss_Mesh);

	// vtkUnstructuredGrid *GaussMesh = vtkUnstructuredGrid::SafeDownCast(outInfo2->Get(vtkDataObject::DATA_OBJECT()));

	// /////////////////////////////////////////////////////////////////////// Gauss Point Visualization ///////////////////////////////////////////////////////////////

	// // Add a quadrature scheme dictionary to the data set. This filter is
	// // solely for our convinience. Typically we would expect that users
	// // provide there own in XML format and use the readers or to generate
	// // them on the fly.

	//   vtkUnstructuredGrid *input=0;
	// vtkSmartPointer<vtkQuadratureSchemeDictionaryGenerator> dictGen = vtkSmartPointer<vtkQuadratureSchemeDictionaryGenerator>::New();
	// dictGen->SetInputData(UGrid);

	// // Interpolate fields to the quadrature points. This generates new field data
	// // arrays, but not a set of points.
	// vtkSmartPointer<vtkQuadraturePointInterpolator> fieldInterp =vtkSmartPointer<vtkQuadraturePointInterpolator>::New();
	// fieldInterp->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "QuadratureOffset");
	// fieldInterp->SetInputConnection(dictGen->GetOutputPort());

	// // Generate the quadrature point set using a specific array as point data.
	// vtkSmartPointer<vtkQuadraturePointsGenerator> pointGen = vtkSmartPointer<vtkQuadraturePointsGenerator>::New();
	// pointGen->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "QuadratureOffset");
	// // pointGen->SetInputConnection(thresholder->GetOutputPort());
	// vtkPolyData *Coutput=vtkPolyData::SafeDownCast(pointGen->GetOutput());
	// pointGen->Update();
	// const char* activeScalars = "pressure";
	// UGrid->GetPointData()->SetActiveScalars(activeScalars);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// cout << "Reached Here" << endl;
	// Build_Gauss_Attributes();
	// Build_Node_Attributes();


	// UGrid_Current_Mesh->ShallowCopy(UGrid_Node_Mesh);
	// GaussMesh->ShallowCopy(UGrid_Gauss_Mesh);
	// output->ShallowCopy(GeneralMesh);

	// cout << "Reached Here" << endl;
	return 1;
}

int pvESSI::RequestInformation( vtkInformation *request, vtkInformationVector **vtkNotUsed(inVec), vtkInformationVector* outVec){

	vtkInformation* Node_Mesh = outVec->GetInformationObject(0);
	vtkInformation* Gauss_Mesh = outVec->GetInformationObject(1);

	//////////////////////////////////////////////////////// Reading Hdf5 File ///////////////////////////////////////////////////////////////////////////////////////

	H5File file = H5File(this->FileName, H5F_ACC_RDWR );
	int Build_Map_Status =1;

	////////////////////////////////////////////////////// Reading General Outside Data /////////////////////////////////////////////////////////////////////////////

	DataSet No_of_TimeSteps_DataSet = file.openDataSet("Number_of_Time_Steps");
	int No_of_TimeSteps[1]; No_of_TimeSteps_DataSet.read( No_of_TimeSteps, PredType::NATIVE_INT);
	this->Number_Of_Time_Steps = No_of_TimeSteps[0];

	DataSet Number_of_Elements_DataSet = file.openDataSet("Number_of_Elements");
	int Number_of_Elements[1]; Number_of_Elements_DataSet.read( Number_of_Elements, PredType::NATIVE_INT);
	this->Number_of_Elements = Number_of_Elements[0];

	DataSet Number_of_Nodes_DataSet = file.openDataSet("Number_of_Nodes");
	int Number_of_Nodes[1]; Number_of_Nodes_DataSet.read( Number_of_Nodes, PredType::NATIVE_INT);
	this->Number_of_Nodes = Number_of_Nodes[0];

	DataSet Whether_Maps_Build_DataSet;
	try { Whether_Maps_Build_DataSet = file.openDataSet("Whether_Maps_Build");}
	catch( FileIException error ){cout << " Dataset is not found." << endl;Build_Map_Status =0; }
	// catch failure caused by the DataSet operations
	catch( DataSetIException error ){cout << "DataSetIException" << endl;}
	// catch failure caused by the DataSpace operations
	catch( DataSpaceIException error ){ cout << "DataSpaceIException" << endl;}
	// catch failure caused by the DataSpace operations
	catch( DataTypeIException error ){ cout << " DataTypeIException" << endl;}

   	file.close(); file = H5File(this->FileName, H5F_ACC_RDONLY ); /// opening in Readonly mode

   	if(!Build_Map_Status){
   		Build_Maps();
   	}
   	else{

	   	DataSet Node_Index_to_Coordinates_DataSet = file.openDataSet("Model/Nodes/Index_to_Coordinates");
		DataSpace Node_Index_to_Coordinates_DataSpace =Node_Index_to_Coordinates_DataSet.getSpace(); hsize_t dims_out[1]; int ndims = Node_Index_to_Coordinates_DataSpace.getSimpleExtentDims( dims_out, NULL);
		this->Pseudo_Number_of_Nodes = dims_out[0];

		DataSet Element_Index_to_Connectivity_DataSet = file.openDataSet("Model/Elements/Index_to_Connectivity");
		DataSpace Element_Index_to_Connectivity_DataSpace =Element_Index_to_Connectivity_DataSet.getSpace(); ndims = Element_Index_to_Connectivity_DataSpace.getSimpleExtentDims( dims_out, NULL);
		this->Pseudo_Number_of_Elements = dims_out[0];
   }

	///////////////////////////////////////////////////////////////// Evaluating Number of Gauss Points ////////////////////////////////////////////////////////////////////////////////////////////

   	DataSet Element_Number_of_Gauss_Points_DataSet = file.openDataSet("Model/Elements/Number_of_Gauss_Points");
	int Element_Number_of_Gauss_Points[Pseudo_Number_of_Elements]; Element_Number_of_Gauss_Points_DataSet.read( Element_Number_of_Gauss_Points, PredType::NATIVE_INT);
	this->Number_of_Gauss_Nodes = 0;

    for (int i=0; i < this->Pseudo_Number_of_Elements; i++){
        if (Element_Number_of_Gauss_Points[i] > 0){
            this->Number_of_Gauss_Nodes += Element_Number_of_Gauss_Points[i];
        }
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	DataSet Time_DataSet = file.openDataSet("time");
	this->Time = new double[Number_Of_Time_Steps]; Time_DataSet.read( Time, PredType::NATIVE_DOUBLE);

	Build_Time_Map(); double Time_range[2]={Time[0],Time[Number_Of_Time_Steps-1]};

	
	Node_Mesh->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),Time, No_of_TimeSteps[0]);
	Node_Mesh->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),Time_range,2);

	// Gauss_Mesh->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),Time, No_of_TimeSteps[0]);
	// Gauss_Mesh->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),Time_range,2);

	// outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),extent,6);
	// outInfo->Set(vtkAlgorithm::CAN_PRODUCE_SUB_EXTENT(),1);


	///////////////////////////////////////////////// Assign Key /////////////////////////////////////////////////////
	// vtkSmartPointer<vtkInformationQuadratureSchemeDefinitionVectorKey> Quadrature_Definition_Vector = vtkSmartPointer<vtkInformationQuadratureSchemeDefinitionVectorKey> ::New();

	// vtkSmartPointer<vtkQuadratureSchemeDefinition> Quadrature_Definition = vtkSmartPointer<vtkQuadratureSchemeDefinition>::New();
	// Quadrature_Definition->Initialize(VTK_HEXAHEDRON,8,8,W_QQ_64_A);

	file.close();

cout << "After Here " << endl;


	return 1;
}
 
void pvESSI::PrintSelf(ostream& os, vtkIndent indent){

	this->Superclass::PrintSelf(os,indent);
	os << indent << "File Name: " << (this->FileName ? this->FileName : "(none)") << "\n";
	return;
}

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

void pvESSI::Build_Node_Attributes(){

	/////////////////////////////////////////////////////////////////// Reading a HDF5 file /////////////////////////////////////////////////////////////////////////////////////////////////

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	DataSet Node_Map_DataSet = file.openDataSet("Maps/NodeMap");
	int *Node_Map; Node_Map = (int*) malloc((Number_of_Nodes) * sizeof(int)); Node_Map_DataSet.read( Node_Map, PredType::NATIVE_INT);

	DataSet Element_Map_DataSet = file.openDataSet("Maps/ElementMap");
	int *Element_Map; Element_Map = (int*) malloc((Number_of_Elements) * sizeof(int)); Element_Map_DataSet.read( Element_Map, PredType::NATIVE_INT);

	DataSet Node_Index_to_Coordinates_DataSet = file.openDataSet("Model/Nodes/Index_to_Coordinates");
	int *Node_Index_to_Coordinates; Node_Index_to_Coordinates = (int*) malloc((Pseudo_Number_of_Nodes) * sizeof(int)); Node_Index_to_Coordinates_DataSet.read( Node_Index_to_Coordinates, PredType::NATIVE_INT);

	DataSet Node_Index_to_Generalized_Displacements_DataSet = file.openDataSet("Model/Nodes/Index_to_Generalized_Displacements");
	int *Node_Index_to_Generalized_Displacements; Node_Index_to_Generalized_Displacements = (int*) malloc((Pseudo_Number_of_Nodes) * sizeof(int)); Node_Index_to_Generalized_Displacements_DataSet.read( Node_Index_to_Generalized_Displacements, PredType::NATIVE_INT);

	///////////////////////////////////////////  Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	
	
    DataSet Node_Generalized_Displacements_DataSet = file.openDataSet("Model/Nodes/Generalized_Displacements"); hsize_t dims_out[2]; int ndims;
	DataSpace Node_Generalized_Displacements_DataSpace = Node_Generalized_Displacements_DataSet.getSpace(); ndims = Node_Generalized_Displacements_DataSpace.getSimpleExtentDims( dims_out, NULL);
	hsize_t col_dims[1];  col_dims[0] = dims_out[0]; DataSpace mspace2( 1, col_dims ); 
	hsize_t offset[2] = { 0, Node_Mesh_Current_Time }; hsize_t  count[2] = { dims_out[0], 1 }; /* Define the column (hyperslab) to read. */
    float *Node_Generalized_Displacements; Node_Generalized_Displacements = (float*) malloc(dims_out[0] * sizeof(float)); // buffer for column to be read and defining hyperslab and read.
    Node_Generalized_Displacements_DataSpace.selectHyperslab( H5S_SELECT_SET, count, offset ); Node_Generalized_Displacements_DataSet.read( Node_Generalized_Displacements, PredType::NATIVE_FLOAT, mspace2, Node_Generalized_Displacements_DataSpace );
 	
 	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	vtkSmartPointer<vtkFloatArray> vtk_Generalized_Displacements = vtkSmartPointer<vtkFloatArray>::New(); 
	vtkSmartPointer<vtkFloatArray> vtk_Generalized_Velocity = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> vtk_Generalized_Acceleration = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> Elastic_Strain_Tensor = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> Plastic_Strain_Tensor = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> Stress_Tensor = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> Material_Properties = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> Total_Energy = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> Incremental_Energy = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkIntArray> Node_Label = vtkSmartPointer<vtkIntArray> ::New();
	vtkSmartPointer<vtkIntArray> Element_Label = vtkSmartPointer<vtkIntArray> ::New();

 	SetMetaArrays( vtk_Generalized_Displacements, vtk_Generalized_Velocity, vtk_Generalized_Acceleration, Elastic_Strain_Tensor, Plastic_Strain_Tensor, Stress_Tensor, Material_Properties, Total_Energy, Incremental_Energy, Node_Label, Element_Label);
 	
 	/////////////////////////////////////////////////////////////// DataSets Visulization at Nodes //////////////////////////////////////////////////////////////////////////////////////
 	int Node_No; 
	for (int i = 0; i < Number_of_Nodes; i++){

			Node_No = Node_Map[i]; Node_Label -> InsertValue(i,Node_No);

			float tuple[3]={Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[Node_No]],Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[Node_No]+1],Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[Node_No]+2]};	

			vtk_Generalized_Displacements->InsertTypedTuple(i,tuple);


			/*************************************************************** Displacement, Acceleration and Velocity -> Calculation Formulae **********************************************************/

			// Displacement(t) = D_t;
			// Acceleration(t) = (1/del.t)^2*{ D_(t+del.t) - 2*D_t + D_(t-del.t)};
			// Acceleration(t) = 1/2/del.t*{ D_(t+del.t) - D_(t-del.t)};

			/***************************************************************** Add Acceleration and Velocity Arrays Here ****************************************/
				
			/****************************************************************************************************************************************************/

	}

	UGrid_Node_Mesh->GetPointData()->AddArray(vtk_Generalized_Displacements);
	UGrid_Node_Mesh->GetPointData()->AddArray(Node_Label);

	// /////////////////////////////////////////////////////////////////////// DataSet Visualization in Cell   ///////////////////////////////////////////////////////////////
 // 	int Element_No; 
	// for (int i = 0; i < Number_of_Elements; i++){

	// 		Element_No = Element_Map[i]; Element_Label -> InsertValue(i,Element_No);
	// }

	// UGrid_Node_Mesh->GetCellData()->AddArray(Element_Label);

	// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	file.close();

  	return;
}

void pvESSI::Build_Gauss_Attributes(){

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	//////////////////////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet Element_Map_DataSet = file.openDataSet("Maps/ElementMap");
	int *Element_Map; Element_Map = (int*) malloc((Number_of_Elements) * sizeof(int)); Element_Map_DataSet.read( Element_Map, PredType::NATIVE_INT);

	DataSet Element_Index_to_Outputs_DataSet = file.openDataSet("Model/Elements/Index_to_Outputs");
	int Element_Index_to_Outputs[Number_of_Elements+1]; Element_Index_to_Outputs_DataSet.read( Element_Index_to_Outputs, PredType::NATIVE_INT);

	DataSet Element_Number_of_Gauss_Points_DataSet = file.openDataSet("Model/Elements/Number_of_Gauss_Points");
	int Element_Number_of_Gauss_Points[Number_of_Elements+1]; Element_Number_of_Gauss_Points_DataSet.read( Element_Number_of_Gauss_Points, PredType::NATIVE_INT);

	DataSet Element_Outputs_DataSet = file.openDataSet("Model/Elements/Outputs"); hsize_t dims_out[2]; int ndims;
	DataSpace Element_Outputs_DataSpace = Element_Outputs_DataSet.getSpace(); ndims = Element_Outputs_DataSpace.getSimpleExtentDims( dims_out, NULL);
	hsize_t col_dims[1];  col_dims[0] = dims_out[0]; DataSpace mspace1( 1, col_dims );
    hsize_t offset[2] = { 0, Gauss_Mesh_Current_Time}; hsize_t  count[2] = { dims_out[0], 1 }; /* Define the column (hyperslab) to read. */
    float *Element_Outputs; Element_Outputs = (float*) malloc(dims_out[0] * sizeof(float)); // buffer for column to be read and defining hyperslab and read.
    Element_Outputs_DataSpace.selectHyperslab( H5S_SELECT_SET, count, offset ); Element_Outputs_DataSet.read( Element_Outputs, PredType::NATIVE_FLOAT, mspace1, Element_Outputs_DataSpace );

	///////////////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	vtkSmartPointer<vtkFloatArray> vtk_Generalized_Displacements = vtkSmartPointer<vtkFloatArray>::New(); 
	vtkSmartPointer<vtkFloatArray> vtk_Generalized_Velocity = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> vtk_Generalized_Acceleration = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> Elastic_Strain_Tensor = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> Plastic_Strain_Tensor = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> Stress_Tensor = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> Material_Properties = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> Total_Energy = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkFloatArray> Incremental_Energy = vtkSmartPointer<vtkFloatArray> ::New();
	vtkSmartPointer<vtkIntArray> Node_Label = vtkSmartPointer<vtkIntArray> ::New();
	vtkSmartPointer<vtkIntArray> Element_Label = vtkSmartPointer<vtkIntArray> ::New();

 	SetMetaArrays( vtk_Generalized_Displacements, vtk_Generalized_Velocity, vtk_Generalized_Acceleration, Elastic_Strain_Tensor, Plastic_Strain_Tensor, Stress_Tensor, Material_Properties, Total_Energy, Incremental_Energy, Node_Label, Element_Label);
 	int Index_of_Gauss_Nodes =0, element_no, No_of_Element_Gauss_Nodes,Output_Index;

	for (int i = 0; i < Number_of_Elements; i++){

		element_no = Element_Map[i];

		No_of_Element_Gauss_Nodes = Element_Number_of_Gauss_Points[element_no];
		Output_Index = Element_Index_to_Outputs[element_no];

		for(int j=0; j< No_of_Element_Gauss_Nodes ; j++){

			float El_Strain_Tuple[9] ={Element_Outputs[Output_Index],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+2],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+5]};
			Output_Index = Output_Index+6;
			float Pl_Strain_Tuple[9] ={Element_Outputs[Output_Index],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+2],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+5]};
			Output_Index = Output_Index+6;
			float Stress_Tuple[9] ={Element_Outputs[Output_Index],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+2],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+5]};
			Output_Index = Output_Index+6;
			Elastic_Strain_Tensor->InsertTypedTuple (Index_of_Gauss_Nodes,El_Strain_Tuple);
			Plastic_Strain_Tensor->InsertTypedTuple (Index_of_Gauss_Nodes,Pl_Strain_Tuple);
			Stress_Tensor->InsertTypedTuple (Index_of_Gauss_Nodes,Stress_Tuple);

		}

		Index_of_Gauss_Nodes+=1;

	}

	UGrid_Gauss_Mesh->GetCellData()->AddArray(Stress_Tensor);
  	UGrid_Gauss_Mesh->GetCellData()->AddArray(Plastic_Strain_Tensor);
  	UGrid_Gauss_Mesh->GetCellData()->AddArray(Elastic_Strain_Tensor);

	// freeing heap allocation for variables 

	// free(Element_Outputs);
	// free(Element_Index_to_Connectivity);
	// free(Element_Index_to_Gauss_Point_Coordinates);
	// free(Element_Index_to_Outputs);
	// free(Element_No_of_Output_Fields);
	// free(Element_Number_of_Gauss_Points);
	// free(Element_Number_of_Nodes);
	// free(Element_Connectivity);
	// free(Element_Gauss_Point_Coordinates);
	// free(Node_Coordinates);
	// free(Node_Index_to_Coordinates);
	// free(Node_Index_to_Generalized_Displacements);
	// free(Node_No_of_DOFs);
	// // free(Element_Outputs);
	// // free(Node_Generalized_Displacements);

	file.close();

	return;
}


void pvESSI::Get_Node_Mesh(){

	/////////////////////////////////////////// Reading a HDF5 file /////////////////////////////////////////////////////////////////////////////////////////////////

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	//////////////////////////////////////////////////// Reading Map Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet Node_Map_DataSet = file.openDataSet("Maps/NodeMap");
	int *Node_Map; Node_Map = (int*) malloc((Number_of_Nodes) * sizeof(int)); Node_Map_DataSet.read( Node_Map, PredType::NATIVE_INT);

	DataSet Inverse_Node_Map_DataSet = file.openDataSet("Maps/InverseNodeMap");
	int *Inverse_Node_Map; Inverse_Node_Map = (int*) malloc((Pseudo_Number_of_Nodes) * sizeof(int)); Inverse_Node_Map_DataSet.read( Inverse_Node_Map, PredType::NATIVE_INT);

	DataSet Element_Map_DataSet = file.openDataSet("Maps/ElementMap");
	int *Element_Map; Element_Map = (int*) malloc((Number_of_Elements) * sizeof(int)); Element_Map_DataSet.read( Element_Map, PredType::NATIVE_INT);

	//////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet Element_Index_to_Connectivity_DataSet = file.openDataSet("Model/Elements/Index_to_Connectivity");
	int *Element_Index_to_Connectivity; Element_Index_to_Connectivity = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int)); Element_Index_to_Connectivity_DataSet.read( Element_Index_to_Connectivity, PredType::NATIVE_INT);

	DataSet Element_Class_Tags_DataSet = file.openDataSet("Model/Elements/Class_Tags");
	int *Element_Class_Tags; Element_Class_Tags = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int)); Element_Class_Tags_DataSet.read( Element_Class_Tags, PredType::NATIVE_INT);

	DataSet Element_Index_to_Gauss_Point_Coordinates_DataSet = file.openDataSet("Model/Elements/Index_to_Gauss_Point_Coordinates");
	int *Element_Index_to_Gauss_Point_Coordinates; Element_Index_to_Gauss_Point_Coordinates = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int)); Element_Index_to_Gauss_Point_Coordinates_DataSet.read( Element_Index_to_Gauss_Point_Coordinates, PredType::NATIVE_INT);

	DataSet Element_Number_of_Gauss_Points_DataSet = file.openDataSet("Model/Elements/Number_of_Gauss_Points");
	int *Element_Number_of_Gauss_Points; Element_Number_of_Gauss_Points = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int)); Element_Number_of_Gauss_Points_DataSet.read( Element_Number_of_Gauss_Points, PredType::NATIVE_INT);

	DataSet Element_Number_of_Nodes_DataSet = file.openDataSet("Model/Elements/Number_of_Nodes");
	int *Element_Number_of_Nodes; Element_Number_of_Nodes = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int)); Element_Number_of_Nodes_DataSet.read( Element_Number_of_Nodes, PredType::NATIVE_INT);

	DataSet Element_Connectivity_DataSet = file.openDataSet("Model/Elements/Connectivity");
	DataSpace Element_Connectivity_DataSpace = Element_Connectivity_DataSet.getSpace(); hsize_t dims_out[1]; int ndims = Element_Connectivity_DataSpace.getSimpleExtentDims( dims_out, NULL);
	int *Element_Connectivity; Element_Connectivity = (int*) malloc(dims_out[0] * sizeof(int));	Element_Connectivity_DataSet.read( Element_Connectivity, PredType::NATIVE_INT);

	DataSet Element_Gauss_Point_Coordinates_DataSet = file.openDataSet("Model/Elements/Gauss_Point_Coordinates");
	DataSpace Element_Gauss_Point_Coordinates_DataSpace = Element_Gauss_Point_Coordinates_DataSet.getSpace(); ndims = Element_Gauss_Point_Coordinates_DataSpace.getSimpleExtentDims( dims_out, NULL);
	float *Element_Gauss_Point_Coordinates; Element_Gauss_Point_Coordinates = (float*) malloc(dims_out[0] * sizeof(float));	Element_Gauss_Point_Coordinates_DataSet.read( Element_Gauss_Point_Coordinates, PredType::NATIVE_FLOAT);

	///////////////////////////////////////////////// Reading Nodes dataset //////////////////////////////////////////////////////////////////////////////////////////////

	DataSet Node_Coordinates_DataSet = file.openDataSet("Model/Nodes/Coordinates");
	float *Node_Coordinates; Node_Coordinates = (float*) malloc((Number_of_Nodes*3) * sizeof(float)); Node_Coordinates_DataSet.read( Node_Coordinates, PredType::NATIVE_FLOAT);

	DataSet Node_Index_to_Coordinates_DataSet = file.openDataSet("Model/Nodes/Index_to_Coordinates");
	int *Node_Index_to_Coordinates; Node_Index_to_Coordinates = (int*) malloc((Pseudo_Number_of_Nodes) * sizeof(int)); Node_Index_to_Coordinates_DataSet.read( Node_Index_to_Coordinates, PredType::NATIVE_INT);

	/////////////////////////////////////////////////// Building Nodes ///////////////////////////////////////////////////////////////////////////////////////////

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	points->SetNumberOfPoints(Number_of_Nodes);
	UGrid_Node_Mesh->Allocate(Number_of_Elements);
	int node_no ;
	for (int i = 0; i < Number_of_Nodes; i++){
		node_no = Node_Map[i];
		points->InsertPoint(i, Node_Coordinates[Node_Index_to_Coordinates[node_no]],Node_Coordinates[Node_Index_to_Coordinates[node_no]+1],Node_Coordinates[Node_Index_to_Coordinates[node_no]+2]);
	}

	UGrid_Node_Mesh->SetPoints(points);

	///////////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	int connectivity_index,No_of_Element_Nodes,Cell_Type,element_no; bool status=true;

	for (int i = 0; i < Number_of_Elements; i++){

		element_no = Element_Map[i];

		connectivity_index = Element_Index_to_Connectivity[element_no];
		No_of_Element_Nodes = Element_Number_of_Nodes[element_no];

		vtkIdType Vertices[No_of_Element_Nodes];
		Cell_Type = ESSI_to_VTK_Element.find(No_of_Element_Nodes)->second;
		std::vector<int> Nodes_Connectivity_Order = ESSI_to_VTK_Connectivity.find(No_of_Element_Nodes)->second;

		for(int j=0; j<No_of_Element_Nodes ; j++){
			Vertices[j] = Inverse_Node_Map[Element_Connectivity[connectivity_index+Nodes_Connectivity_Order[j]]];
		}

		UGrid_Node_Mesh->InsertNextCell(Cell_Type, No_of_Element_Nodes, Vertices);

	}

	Whether_Node_Mesh_Build=1;

	file.close();

	return;
	
}

void pvESSI::Get_Gauss_Mesh(){

	/////////////////////////////////////////// Reading a HDF5 file /////////////////////////////////////////////////////////////////////////////////////////////////

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	/////////////////////////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet Element_Map_DataSet = file.openDataSet("Maps/ElementMap");
	int *Element_Map; Element_Map = (int*) malloc((Number_of_Elements) * sizeof(int)); Element_Map_DataSet.read( Element_Map, PredType::NATIVE_INT);

	DataSet Element_Index_to_Gauss_Point_Coordinates_DataSet = file.openDataSet("Model/Elements/Index_to_Gauss_Point_Coordinates");
	int *Element_Index_to_Gauss_Point_Coordinates; Element_Index_to_Gauss_Point_Coordinates = (int*) malloc((Number_of_Elements+1) * sizeof(int)); Element_Index_to_Gauss_Point_Coordinates_DataSet.read( Element_Index_to_Gauss_Point_Coordinates, PredType::NATIVE_INT);

	DataSet Element_Number_of_Gauss_Points_DataSet = file.openDataSet("Model/Elements/Number_of_Gauss_Points");
	int *Element_Number_of_Gauss_Points; Element_Number_of_Gauss_Points = (int*) malloc((Number_of_Elements+1) * sizeof(int)); Element_Number_of_Gauss_Points_DataSet.read( Element_Number_of_Gauss_Points, PredType::NATIVE_INT);

 	DataSet Element_Gauss_Point_Coordinates_DataSet = file.openDataSet("Model/Elements/Gauss_Point_Coordinates");
	DataSpace Element_Gauss_Point_Coordinates_DataSpace = Element_Gauss_Point_Coordinates_DataSet.getSpace(); hsize_t dims_out[1]; int ndims = Element_Gauss_Point_Coordinates_DataSpace.getSimpleExtentDims( dims_out, NULL);
	float *Element_Gauss_Point_Coordinates; Element_Gauss_Point_Coordinates = (float*) malloc(dims_out[0] * sizeof(float));	Element_Gauss_Point_Coordinates_DataSet.read( Element_Gauss_Point_Coordinates, PredType::NATIVE_FLOAT);

	////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(Number_of_Gauss_Nodes); int Index_of_Gauss_Nodes=0, element_no;
	int No_of_Element_Gauss_Nodes;

	for (int i=0; i < Number_of_Elements; i++){
		element_no = Element_Map[i];
		No_of_Element_Gauss_Nodes = Element_Number_of_Gauss_Points[element_no];
		for(int j=0; j<No_of_Element_Gauss_Nodes ; j++){
			points->InsertPoint(Index_of_Gauss_Nodes, Element_Gauss_Point_Coordinates[Element_Index_to_Gauss_Point_Coordinates[element_no]+3*j],Element_Gauss_Point_Coordinates[Element_Index_to_Gauss_Point_Coordinates[element_no]+3*j+1],Element_Gauss_Point_Coordinates[Element_Index_to_Gauss_Point_Coordinates[element_no]+3*j+2]);
			Index_of_Gauss_Nodes = Index_of_Gauss_Nodes+1;
		}

	}

	UGrid_Gauss_Mesh->SetPoints(points);
	Whether_Gauss_Mesh_Build=1;

	file.close();

	return;
	
}

void pvESSI::SetMetaArrays( vtkSmartPointer<vtkFloatArray> &vtk_Generalized_Displacements, vtkSmartPointer<vtkFloatArray> &vtk_Generalized_Velocity, vtkSmartPointer<vtkFloatArray> &vtk_Generalized_Acceleration, vtkSmartPointer<vtkFloatArray> &Elastic_Strain_Tensor, vtkSmartPointer<vtkFloatArray> &Plastic_Strain_Tensor, vtkSmartPointer<vtkFloatArray> &Stress_Tensor, vtkSmartPointer<vtkFloatArray> &Material_Properties, vtkSmartPointer<vtkFloatArray> &Total_Energy, vtkSmartPointer<vtkFloatArray> &Incremental_Energy, vtkSmartPointer<vtkIntArray> &Node_Label, vtkSmartPointer<vtkIntArray> &Element_Label ){

	/*************************************** Mesh Meta Data ************************************************************/

	vtk_Generalized_Displacements->SetName("Generalized_Displacements");
	vtk_Generalized_Displacements->SetNumberOfComponents(3);
	vtk_Generalized_Displacements->SetComponentName(0,"X-axis");
	vtk_Generalized_Displacements->SetComponentName(1,"Y-axis");
	vtk_Generalized_Displacements->SetComponentName(2,"Z-axis");

	vtk_Generalized_Velocity->SetName("Velocity");
	vtk_Generalized_Velocity->SetNumberOfComponents(3);
	vtk_Generalized_Velocity->SetComponentName(0,"X-axis");
	vtk_Generalized_Velocity->SetComponentName(1,"Y-axis");
	vtk_Generalized_Velocity->SetComponentName(2,"Z-axis");

	vtk_Generalized_Acceleration->SetName("Acceleration");
	vtk_Generalized_Acceleration->SetNumberOfComponents(3);
	vtk_Generalized_Acceleration->SetComponentName(0,"X-axis");
	vtk_Generalized_Acceleration->SetComponentName(1,"Y-axis");
	vtk_Generalized_Acceleration->SetComponentName(2,"Z-axis");

	Elastic_Strain_Tensor->SetName("Elastic_Strain");
	Elastic_Strain_Tensor->SetNumberOfComponents(9);
	Elastic_Strain_Tensor->SetComponentName(0,"El-Strain-xx");
	Elastic_Strain_Tensor->SetComponentName(1,"El-Strain-xy");
	Elastic_Strain_Tensor->SetComponentName(2,"El-Strain-xz");
	Elastic_Strain_Tensor->SetComponentName(3,"El-Strain-yx");
	Elastic_Strain_Tensor->SetComponentName(4,"El-Strain-yy");
	Elastic_Strain_Tensor->SetComponentName(5,"El-Strain-yz");
	Elastic_Strain_Tensor->SetComponentName(6,"El-Strain-zx");
	Elastic_Strain_Tensor->SetComponentName(7,"El-Strain-zy");
	Elastic_Strain_Tensor->SetComponentName(8,"El-Strain-zz");

	Plastic_Strain_Tensor->SetName("Plastic_Strain");
	Plastic_Strain_Tensor->SetNumberOfComponents(9);
	Plastic_Strain_Tensor->SetComponentName(0,"Pl-Strain-xx");
	Plastic_Strain_Tensor->SetComponentName(1,"Pl-Strain-xy");
	Plastic_Strain_Tensor->SetComponentName(2,"Pl-Strain-xz");
	Plastic_Strain_Tensor->SetComponentName(3,"Pl-Strain-yx");
	Plastic_Strain_Tensor->SetComponentName(4,"Pl-Strain-yy");
	Plastic_Strain_Tensor->SetComponentName(5,"Pl-Strain-yz");
	Plastic_Strain_Tensor->SetComponentName(6,"Pl-Strain-zx");
	Plastic_Strain_Tensor->SetComponentName(7,"Pl-Strain-zy");
	Plastic_Strain_Tensor->SetComponentName(8,"Pl-Strain-zz");	

	Stress_Tensor->SetName("Stress");
	Stress_Tensor->SetNumberOfComponents(9);
	Stress_Tensor->SetComponentName(0,"Stress-xx");
	Stress_Tensor->SetComponentName(1,"Stress-xy");
	Stress_Tensor->SetComponentName(2,"Stress-xz");
	Stress_Tensor->SetComponentName(3,"Stress-yx");
	Stress_Tensor->SetComponentName(4,"Stress-yy");
	Stress_Tensor->SetComponentName(5,"Stress-yz");
	Stress_Tensor->SetComponentName(6,"Stress-zx");
	Stress_Tensor->SetComponentName(7,"Stress-zy");
	Stress_Tensor->SetComponentName(8,"Stress-zz");	

	Material_Properties->SetName("Material_Properties");
	Material_Properties->SetNumberOfComponents(1);

	Total_Energy->SetName("Total-Energy");
	Total_Energy->SetNumberOfComponents(1);

	Incremental_Energy->SetName("Incremental_Energy");
	Incremental_Energy->SetNumberOfComponents(1);

	Node_Label->SetName("Node_No");
	Node_Label->SetNumberOfComponents(1);

	Element_Label->SetName("Element_No");
	Element_Label->SetNumberOfComponents(1);

	return;

	/*************************************************************************************************************************************************/
}

void pvESSI::initialize(){

	Node_Mesh_Current_Time=0;
	Gauss_Mesh_Current_Time=0;
  	Node_Mesh_Previous_Time=-1;
  	Gauss_Mesh_Previous_Time=-1; 
	Display_Node_Mesh=1;
	Display_Gauss_Mesh=0;
	Display_All_Mesh=0;
	Whether_Node_Mesh_Build=0;
	Whether_Gauss_Mesh_Build=0;
	Whether_All_Mesh_Build=0;
	// this->DebugOn();
	// cout << "xascaS" << endl;
}

void pvESSI::Build_Maps(){

	H5File file = H5File(this->FileName, H5F_ACC_RDWR );

	DataSet Node_Index_to_Coordinates_DataSet = file.openDataSet("Model/Nodes/Index_to_Coordinates");
	DataSpace Node_Index_to_Coordinates_DataSpace =Node_Index_to_Coordinates_DataSet.getSpace(); hsize_t dims_out[1]; int ndims = Node_Index_to_Coordinates_DataSpace.getSimpleExtentDims( dims_out, NULL);
	this->Pseudo_Number_of_Nodes = dims_out[0];
	int *Node_Index_to_Coordinates; Node_Index_to_Coordinates = (int*) malloc((Pseudo_Number_of_Nodes) * sizeof(int)); Node_Index_to_Coordinates_DataSet.read( Node_Index_to_Coordinates, PredType::NATIVE_INT);

	DataSet Element_Index_to_Connectivity_DataSet = file.openDataSet("Model/Elements/Index_to_Connectivity");
	DataSpace Element_Index_to_Connectivity_DataSpace =Element_Index_to_Connectivity_DataSet.getSpace(); ndims = Element_Index_to_Connectivity_DataSpace.getSimpleExtentDims( dims_out, NULL);
	this->Pseudo_Number_of_Elements = dims_out[0];
	int *Element_Index_to_Connectivity; Element_Index_to_Connectivity = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int)); Element_Index_to_Connectivity_DataSet.read( Element_Index_to_Connectivity, PredType::NATIVE_INT);

	int *Node_Map; Node_Map = (int*) malloc((Number_of_Nodes) * sizeof(int));
	int *Element_Map; Element_Map = (int*) malloc((Number_of_Elements) * sizeof(int));
	int *Inverse_Node_Map; Inverse_Node_Map = (int*) malloc((Pseudo_Number_of_Nodes) * sizeof(int));
	int *Inverse_Element_Map; Inverse_Element_Map = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int));

	int index =0;
	for (int i = 0; i < Pseudo_Number_of_Nodes; i++)
	{
		if(Node_Index_to_Coordinates[i] !=-1){
			Node_Map[index] = i;
			Inverse_Node_Map[i] = index;
			index = index+1;
		}
		else{
			Inverse_Node_Map[i] = -1;
		}
	
	}

	cout << "Number of Nodes " << index << " " << Number_of_Nodes << endl;

	index =0;
	for (int i = 0; i < Pseudo_Number_of_Elements; ++i)
	{
		if(Element_Index_to_Connectivity[i] !=-1){
			Element_Map[index] = i;
			Inverse_Element_Map[i] = index;
			index = index+1;
		}
		else{
			Inverse_Element_Map[i] = -1;
		}
	}

	cout << "Number of Elements " << index << " " << Number_of_Elements << endl;;

try{
	Group group(file.createGroup("Maps"));
	dims_out[0]= Number_of_Nodes; DataSpace Node_Map_DataSpace (1,dims_out);
	DataSet Node_Map_DataSet = file.createDataSet("Maps/NodeMap",PredType::NATIVE_INT,Node_Map_DataSpace);
	Node_Map_DataSet.write(Node_Map, PredType::NATIVE_INT );

	dims_out[0] = Number_of_Elements; DataSpace Element_Map_DataSpace (1,dims_out);
	DataSet Element_Map_DataSet = file.createDataSet("Maps/ElementMap",PredType::NATIVE_INT,Element_Map_DataSpace);
	Element_Map_DataSet.write(Element_Map, PredType::NATIVE_INT );

	dims_out[0] = Pseudo_Number_of_Nodes; DataSpace Inverse_Node_Map_DataSpace (1,dims_out);
	DataSet Inverse_Node_Map_DataSet = file.createDataSet("Maps/InverseNodeMap",PredType::NATIVE_INT,Inverse_Node_Map_DataSpace);
	Inverse_Node_Map_DataSet.write(Inverse_Node_Map,PredType::NATIVE_INT);

	dims_out[0] = Pseudo_Number_of_Elements; DataSpace Inverse_Element_Map_DataSpace (1,dims_out);
	DataSet Inverse_Element_Map_DataSet = file.createDataSet("Maps/InverseElementMap",PredType::NATIVE_INT,Inverse_Element_Map_DataSpace);
	Inverse_Element_Map_DataSet.write(Inverse_Node_Map,PredType::NATIVE_INT);

	dims_out[0] = 1; DataSpace Whether_Maps_Build_DataSpace (1,dims_out);int Whether_Maps_Build[1]; Whether_Maps_Build [0]= 1;
	DataSet Whether_Maps_Build_DataSet = file.createDataSet("Whether_Maps_Build",PredType::NATIVE_INT,Whether_Maps_Build_DataSpace);
	Whether_Maps_Build_DataSet.write(Whether_Maps_Build,PredType::NATIVE_INT);


}

   catch( FileIException error )
   {
      cout << " Dataset is not found." << endl;
   }
   // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
      cout << "DataSetIException" << endl;
   }
   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
      cout << "DataSpaceIException" << endl;
   }
   // catch failure caused by the DataSpace operations
   catch( DataTypeIException error )
   {
     cout << " DataTypeIException" << endl;
   }

	file.close();
}

void pvESSI::Build_Time_Map(){

	for (int i = 0;i<Number_Of_Time_Steps ; i++)
		Time_Map[Time[i]] = i;

	return;
}



