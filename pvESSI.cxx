#include "pvESSI.h"
#include <stdlib.h>

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
// cmake .. -DParaView_DIR=/home/sumeet/Softwares/Apllications/Paraview -DGHOST_BUILD_CDAWEB=OFF
/************************************************************************************************************************************************/

using namespace::H5;

vtkStandardNewMacro(pvESSI);
 
pvESSI::pvESSI(){ 

	this->FileName = NULL;
	this->SetNumberOfInputPorts(0);
	this->SetNumberOfOutputPorts(1);
	Energy_Database_Status=-1;
	mesh_build =0;
	UGrid_Mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
	this->set_VTK_To_ESSI_Elements_Connectivity();
}

int pvESSI::RequestData(vtkInformation *vtkNotUsed(request),vtkInformationVector **vtkNotUsed(inputVector),	vtkInformationVector *outputVector){
 
	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	// outInfo->Print(std::cout);

	if (!mesh_build) this->GetMesh(); 
	
	// if(Energy_Database_Status==-1)
	// 	this->Make_Energy_Database();

	// int extent[6] = {0,-1,0,-1,0,-1};
	// outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);

  	int Ctime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
  	// int Clength = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  	// double* Csteps = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

 	// cout<< "Ctime " << "  " << Ctime << endl;
	// cout<< "Clength " << "  " << Clength << endl;
	// for (int i =0 ; i< Clength ; i++){
	// 	cout << Csteps[i] << "  " << endl;
	// }

	// if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
	//          std::cout << "Update Time" << endl;
	//     if (!outInfo->Has(vtkDataObject::DATA_TIME_STEP()))
	//     		std::cout << "DATA_TIME_STEP" << endl;

	// return outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())
	// int piece, numPieces;
	// piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
	// numPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

	/////////////////////////////// Read in Paralll //////////////////////////////
	// outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
	///////////////////////////////////////////////////////////////////////////////

	// get the ouptut pointer to paraview 
	vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	///////////////////////////////////////////// Reading a HDF5 file /////////////////////////////////////////////////////////////////////////////////////////////////

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	////////////////////////////////////// Reading General Outside Data ////////////////////////////////////////////////////////////////////////////////////////////

	DataSet No_of_Elements_DataSet = file.openDataSet("Number_of_Elements");
	int Number_of_Elements[1]; No_of_Elements_DataSet.read( Number_of_Elements, PredType::NATIVE_INT);

	DataSet Number_of_Nodes_DataSet = file.openDataSet("Number_of_Nodes");
	int Number_of_Nodes[1]; Number_of_Nodes_DataSet.read( Number_of_Nodes, PredType::NATIVE_INT);

	DataSet No_of_TimeSteps_DataSet = file.openDataSet("Number_of_Time_Steps");
	int No_of_TimeSteps[1]; No_of_TimeSteps_DataSet.read( No_of_TimeSteps, PredType::NATIVE_INT);

	//////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet Element_Index_to_Connectivity_DataSet = file.openDataSet("Model/Elements/Index_to_Connectivity");
	int *Element_Index_to_Connectivity; Element_Index_to_Connectivity = (int*) malloc((Number_of_Elements[0]+1) * sizeof(int)); Element_Index_to_Connectivity_DataSet.read( Element_Index_to_Connectivity, PredType::NATIVE_INT);

	DataSet Element_Index_to_Gauss_Point_Coordinates_DataSet = file.openDataSet("Model/Elements/Index_to_Gauss_Point_Coordinates");
	int *Element_Index_to_Gauss_Point_Coordinates; Element_Index_to_Gauss_Point_Coordinates = (int*) malloc((Number_of_Elements[0]+1) * sizeof(int)); Element_Index_to_Gauss_Point_Coordinates_DataSet.read( Element_Index_to_Gauss_Point_Coordinates, PredType::NATIVE_INT);

	DataSet Element_Index_to_Outputs_DataSet = file.openDataSet("Model/Elements/Index_to_Outputs");
	int *Element_Index_to_Outputs; Element_Index_to_Outputs = (int*) malloc((Number_of_Elements[0]+1) * sizeof(int)); Element_Index_to_Outputs_DataSet.read( Element_Index_to_Outputs, PredType::NATIVE_INT);

	DataSet Element_No_of_Output_Fields_DataSet = file.openDataSet("Model/Elements/Number_of_Output_Fields");
	int *Element_No_of_Output_Fields; Element_No_of_Output_Fields = (int*) malloc((Number_of_Elements[0]+1) * sizeof(int)); Element_No_of_Output_Fields_DataSet.read( Element_No_of_Output_Fields, PredType::NATIVE_INT);

	DataSet Element_Number_of_Gauss_Points_DataSet = file.openDataSet("Model/Elements/Number_of_Gauss_Points");
	int *Element_Number_of_Gauss_Points; Element_Number_of_Gauss_Points = (int*) malloc((Number_of_Elements[0]+1) * sizeof(int)); Element_Number_of_Gauss_Points_DataSet.read( Element_Number_of_Gauss_Points, PredType::NATIVE_INT);

	DataSet Element_Number_of_Nodes_DataSet = file.openDataSet("Model/Elements/Number_of_Nodes");
	int *Element_Number_of_Nodes; Element_Number_of_Nodes = (int*) malloc((Number_of_Elements[0]+1) * sizeof(int)); Element_Number_of_Nodes_DataSet.read( Element_Number_of_Nodes, PredType::NATIVE_INT);

	DataSet Element_Connectivity_DataSet = file.openDataSet("Model/Elements/Connectivity");
	DataSpace Element_Connectivity_DataSpace = Element_Connectivity_DataSet.getSpace(); hsize_t dims_out[2]; int ndims = Element_Connectivity_DataSpace.getSimpleExtentDims( dims_out, NULL);
	int *Element_Connectivity; Element_Connectivity = (int*) malloc(dims_out[0] * sizeof(int));	Element_Connectivity_DataSet.read( Element_Connectivity, PredType::NATIVE_INT);

	DataSet Element_Gauss_Point_Coordinates_DataSet = file.openDataSet("Model/Elements/Gauss_Point_Coordinates");
	DataSpace Element_Gauss_Point_Coordinates_DataSpace = Element_Gauss_Point_Coordinates_DataSet.getSpace(); ndims = Element_Gauss_Point_Coordinates_DataSpace.getSimpleExtentDims( dims_out, NULL);
	float *Element_Gauss_Point_Coordinates; Element_Gauss_Point_Coordinates = (float*) malloc(dims_out[0] * sizeof(float));	Element_Gauss_Point_Coordinates_DataSet.read( Element_Gauss_Point_Coordinates, PredType::NATIVE_FLOAT);

	/////////////////////////////////////////////////// Reading Nodes datat //////////////////////////////////////////////////////////////////////////////////////////////

	DataSet Node_Coordinates_DataSet = file.openDataSet("Model/Nodes/Coordinates");
	float *Node_Coordinates; Node_Coordinates = (float*) malloc((Number_of_Nodes[0]*3) * sizeof(float)); Node_Coordinates_DataSet.read( Node_Coordinates, PredType::NATIVE_FLOAT);

	DataSet Node_Index_to_Coordinates_DataSet = file.openDataSet("Model/Nodes/Index_to_Coordinates");
	int *Node_Index_to_Coordinates; Node_Index_to_Coordinates = (int*) malloc((Number_of_Nodes[0]+1) * sizeof(int)); Node_Index_to_Coordinates_DataSet.read( Node_Index_to_Coordinates, PredType::NATIVE_INT);

	DataSet Node_Index_to_Generalized_Displacements_DataSet = file.openDataSet("Model/Nodes/Index_to_Generalized_Displacements");
	int *Node_Index_to_Generalized_Displacements; Node_Index_to_Generalized_Displacements = (int*) malloc((Number_of_Nodes[0]+1) * sizeof(int)); Node_Index_to_Generalized_Displacements_DataSet.read( Node_Index_to_Generalized_Displacements, PredType::NATIVE_INT);

	DataSet Node_No_of_DOFs_DataSet = file.openDataSet("Model/Nodes/Number_of_DOFs");
	int *Node_No_of_DOFs; Node_No_of_DOFs = (int*) malloc((Number_of_Nodes[0]+1) * sizeof(int)); Node_No_of_DOFs_DataSet.read( Node_No_of_DOFs, PredType::NATIVE_INT);

	///////////////////////////////////////////  Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	
	
	DataSet Element_Outputs_DataSet = file.openDataSet("Model/Elements/Outputs");
	DataSpace Element_Outputs_DataSpace = Element_Outputs_DataSet.getSpace(); ndims = Element_Outputs_DataSpace.getSimpleExtentDims( dims_out, NULL);
	hsize_t col_dims[1];  col_dims[0] = dims_out[0]; DataSpace mspace1( 1, col_dims );
    hsize_t offset[2] = { 0, Ctime }; hsize_t  count[2] = { dims_out[0], 1 }; /* Define the column (hyperslab) to read. */
    float *Element_Outputs; Element_Outputs = (float*) malloc(dims_out[0] * sizeof(float)); // buffer for column to be read and defining hyperslab and read.
    Element_Outputs_DataSpace.selectHyperslab( H5S_SELECT_SET, count, offset ); Element_Outputs_DataSet.read( Element_Outputs, PredType::NATIVE_FLOAT, mspace1, Element_Outputs_DataSpace );

    DataSet Node_Generalized_Displacements_DataSet = file.openDataSet("Model/Nodes/Generalized_Displacements");
	DataSpace Node_Generalized_Displacements_DataSpace = Node_Generalized_Displacements_DataSet.getSpace(); ndims = Node_Generalized_Displacements_DataSpace.getSimpleExtentDims( dims_out, NULL);
	col_dims[0] = dims_out[0]; DataSpace mspace2( 1, col_dims ); count[0] = dims_out[0]; /* Define the column (hyperslab) to read. */
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

 	SetMetaArrays( vtk_Generalized_Displacements, vtk_Generalized_Velocity, vtk_Generalized_Acceleration, Elastic_Strain_Tensor, Plastic_Strain_Tensor, Stress_Tensor, Material_Properties, Total_Energy, Incremental_Energy);
	// float *pts  = (float *) points->GetVoidPointer(0);

	for (int i = 0; i <= Number_of_Nodes[0]; i++){

		if (Node_Index_to_Coordinates[i]!=-1){

			float tuple[3]={Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[i]],Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[i]+1],Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[i]+2]};	

			/*************************************************************** Displacement, Acceleration and Velocity -> Calculation Formulae **********************************************************/

			// Displacement(t) = D_t;
			// Acceleration(t) = (1/del.t)^2*{ D_(t+del.t) - 2*D_t + D_(t-del.t)};
			// Acceleration(t) = 1/2/del.t*{ D_(t+del.t) - D_(t-del.t)};

			/***************************************************************** Add Acceleration and Velocity Arrays Here ****************************************/

				vtk_Generalized_Displacements->InsertTupleValue(i,tuple);
				// vtk_Generalized_Velocity->InsertTupleValue(i,tuple);
				// vtk_Generalized_Acceleration->InsertTupleValue(i,tuple);

			/****************************************************************************************************************************************************/
		}
	}

	/////////////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	for (int i = 0; i <= Number_of_Elements[0]; i++){

		int connectivity_index = Element_Index_to_Connectivity[i];

		if(connectivity_index!=-1){
	
			int No_of_Element_Nodes = Element_Number_of_Nodes[i];
			vtkIdType Vertices[No_of_Element_Nodes];

			int Cell_Type = ESSI_to_VTK_Element.find(No_of_Element_Nodes)->second;
			std::vector<int> Nodes_Connectivity_Order = ESSI_to_VTK_Connectivity.find(No_of_Element_Nodes)->second;

			/************************************************** Gauss Calculations **********************************************/

			// int NOF_Gauss_Points = Element_Number_of_Gauss_Points[i];
			// int Output_Index = Element_Index_to_Outputs[i];
			// int NOF_Output_Fields = Element_No_of_Output_Fields[i];
			// float Energy_Increment=0;

			// for(int j=0; j< NOF_Gauss_Points ; j++){
			// 	float El_Strain_Tuple[9] ={Element_Outputs[Output_Index][Ctime],Element_Outputs[Output_Index+1][Ctime],Element_Outputs[Output_Index+3][Ctime],Element_Outputs[Output_Index+1][Ctime],Element_Outputs[Output_Index+2][Ctime],Element_Outputs[Output_Index+4][Ctime],Element_Outputs[Output_Index+3][Ctime],Element_Outputs[Output_Index+4][Ctime],Element_Outputs[Output_Index+5][Ctime]};
			// 	Output_Index = Output_Index+6;
			// 	float Pl_Strain_Tuple[9] ={Element_Outputs[Output_Index][Ctime],Element_Outputs[Output_Index+1][Ctime],Element_Outputs[Output_Index+3][Ctime],Element_Outputs[Output_Index+1][Ctime],Element_Outputs[Output_Index+2][Ctime],Element_Outputs[Output_Index+4][Ctime],Element_Outputs[Output_Index+3][Ctime],Element_Outputs[Output_Index+4][Ctime],Element_Outputs[Output_Index+5][Ctime]};
			// 	Output_Index = Output_Index+6;
			// 	float Stress_Tuple[9] ={Element_Outputs[Output_Index][Ctime],Element_Outputs[Output_Index+1][Ctime],Element_Outputs[Output_Index+3][Ctime],Element_Outputs[Output_Index+1][Ctime],Element_Outputs[Output_Index+2][Ctime],Element_Outputs[Output_Index+4][Ctime],Element_Outputs[Output_Index+3][Ctime],Element_Outputs[Output_Index+4][Ctime],Element_Outputs[Output_Index+5][Ctime]};
			// 	Output_Index = Output_Index+6;
			// 	Elastic_Strain_Tensor->InsertNextTuple(El_Strain_Tuple);
			// 	Plastic_Strain_Tensor->InsertNextTuple(Pl_Strain_Tuple);
			// 	Stress_Tensor->InsertNextTuple(Stress_Tuple);
			// }

			// Incremental_Energy->InsertValue(i,this->Incremental_Energy_Database[Ctime*(Number_of_Elements[0]+1)+i]);
			// Total_Energy->InsertValue(i,this->Total_Energy_Database[Ctime*(Number_of_Elements[0]+1)+i]);

			/********************************************************************************************************************/

		}
	}
	
	/////////////////////////////////////////////////////////////////////// Enable the Data Array  to show in Paraview  ///////////////////////////////////////////////////////////////

	UGrid_Mesh->GetPointData()->AddArray(vtk_Generalized_Displacements);
  // 	UGrid_Mesh->GetPointData()->AddArray(vtk_Generalized_Velocity);
  // 	UGrid_Mesh->GetPointData()->AddArray(vtk_Generalized_Acceleration);
  // 	UGrid_Mesh->GetFieldData()->AddArray(Stress_Tensor);
  // 	UGrid_Mesh->GetFieldData()->AddArray(Plastic_Strain_Tensor);
  // 	UGrid_Mesh->GetFieldData()->AddArray(Elastic_Strain_Tensor);
  // 	UGrid_Mesh->GetCellData()->AddArray(Material_Properties);
  // 	UGrid_Mesh->GetCellData()->AddArray(Total_Energy);
  // 	UGrid_Mesh->GetCellData()->AddArray(Incremental_Energy);

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

	// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	output->DeepCopy(UGrid_Mesh);

	// freeing heap allocation for variables 

	free(Element_Outputs);
	free(Element_Index_to_Connectivity);
	free(Element_Index_to_Gauss_Point_Coordinates);
	free(Element_Index_to_Outputs);
	free(Element_No_of_Output_Fields);
	free(Element_Number_of_Gauss_Points);
	free(Element_Number_of_Nodes);
	free(Element_Connectivity);
	free(Element_Gauss_Point_Coordinates);
	free(Node_Coordinates);
	free(Node_Index_to_Coordinates);
	free(Node_Index_to_Generalized_Displacements);
	free(Node_No_of_DOFs);
	// free(Element_Outputs);
	// free(Node_Generalized_Displacements);

	return 1;
}

int pvESSI::RequestInformation( vtkInformation *request, vtkInformationVector **vtkNotUsed(inVec), vtkInformationVector* outVec){


	vtkInformation* outInfo = outVec->GetInformationObject(0);

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	//////////////////////////////////////// Reading General Outside Data ////////////////////////////////////////////////////////////////////////////////////////////

	DataSet No_of_TimeSteps_DataSet = file.openDataSet("Number_of_Time_Steps");
	int No_of_TimeSteps[1]; No_of_TimeSteps_DataSet.read( No_of_TimeSteps, PredType::NATIVE_INT);

	// DataSet Time_DataSet = file.openDataSet("time");
	// double Time[No_of_TimeSteps[0]]; Time_DataSet.read( Time, PredType::NATIVE_INT);

	this->Number_Of_Time_Steps = No_of_TimeSteps[0];
	double Time_Steps[No_of_TimeSteps[0]]; Time_Steps[0]=0;

	for (int i =0;i<No_of_TimeSteps[0] ;i++)
		Time_Steps[i] = i;

	double Time_range[2]={0,No_of_TimeSteps[0]};

	// // int extent[6]={0,100,0,100,0,100}; 


	// // outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),extent,6);
	outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),Time_Steps, No_of_TimeSteps[0]);
	outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),Time_range,2);
	// // outInfo->Set(vtkAlgorithm::CAN_PRODUCE_SUB_EXTENT(),1);


	///////////////////////////////////////////////// Assign Key /////////////////////////////////////////////////////
	// vtkSmartPointer<vtkInformationQuadratureSchemeDefinitionVectorKey> Quadrature_Definition_Vector = vtkSmartPointer<vtkInformationQuadratureSchemeDefinitionVectorKey> ::New();

	// vtkSmartPointer<vtkQuadratureSchemeDefinition> Quadrature_Definition = vtkSmartPointer<vtkQuadratureSchemeDefinition>::New();
	// Quadrature_Definition->Initialize(VTK_HEXAHEDRON,8,8,W_QQ_64_A);

	return 1;
}
 
void pvESSI::PrintSelf(ostream& os, vtkIndent indent){

	this->Superclass::PrintSelf(os,indent);
	os << indent << "File Name: " << (this->FileName ? this->FileName : "(none)") << "\n";
}

void pvESSI::set_VTK_To_ESSI_Elements_Connectivity(){

	std::vector<int> connectivity_vector;
	// this->TimestepValues[0]=1;

	ESSI_to_VTK_Element[1] = VTK_VERTEX;		ESSI_to_VTK_Element[2]  = VTK_LINE;		ESSI_to_VTK_Element[4] = VTK_QUAD;
	ESSI_to_VTK_Element[8] = VTK_HEXAHEDRON;	/*ESSI_to_VTK_Element[20] = VTK_VERTEX;*/	ESSI_to_VTK_Element[27] = VTK_TRIQUADRATIC_HEXAHEDRON;

	int Node[1] = {0};						connectivity_vector.assign(Node,Node+1);				ESSI_to_VTK_Connectivity[1] = connectivity_vector;	connectivity_vector.clear();
	int Line[2] = {0,1}; 					connectivity_vector.assign(Line,Line+2); 				ESSI_to_VTK_Connectivity[2] = connectivity_vector;	connectivity_vector.clear();
	int Quadrangle[4] = {0,1,2,3}; 			connectivity_vector.assign(Quadrangle,Quadrangle+4); 	ESSI_to_VTK_Connectivity[4] = connectivity_vector;	connectivity_vector.clear();
	int Hexahedron[8] = {4,5,6,7,0,1,2,3}; 	connectivity_vector.assign(Hexahedron,Hexahedron+8);	ESSI_to_VTK_Connectivity[8] = connectivity_vector;	connectivity_vector.clear();
	int TriQudratic_hexahedron[27] = {6,5,4,7,2,1,0,3,13,12,15,14,9,8,11,10,18,17,16,19,23,21,22,24,26,25,20}; 											connectivity_vector.assign(TriQudratic_hexahedron,TriQudratic_hexahedron+27);					ESSI_to_VTK_Connectivity[27]= connectivity_vector;	connectivity_vector.clear();
}

void pvESSI::Make_Energy_Database(){

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	//////////////////////////////////////// Reading General Outside Data ////////////////////////////////////////////////////////////////////////////////////////////

	DataSet No_of_Elements_DataSet = file.openDataSet("Number_of_Elements");
	int Number_of_Elements[1]; No_of_Elements_DataSet.read( Number_of_Elements, PredType::NATIVE_INT);

	DataSet No_of_TimeSteps_DataSet = file.openDataSet("Number_of_Time_Steps");
	int No_of_TimeSteps[1]; No_of_TimeSteps_DataSet.read( No_of_TimeSteps, PredType::NATIVE_INT);

	//////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet Element_Index_to_Connectivity_DataSet = file.openDataSet("Model/Elements/Index_to_Connectivity");
	int Element_Index_to_Connectivity[Number_of_Elements[0]+1]; Element_Index_to_Connectivity_DataSet.read( Element_Index_to_Connectivity, PredType::NATIVE_INT);

	DataSet Element_Index_to_Outputs_DataSet = file.openDataSet("Model/Elements/Index_to_Outputs");
	int Element_Index_to_Outputs[Number_of_Elements[0]+1]; Element_Index_to_Outputs_DataSet.read( Element_Index_to_Outputs, PredType::NATIVE_INT);

	DataSet Element_No_of_Output_Fields_DataSet = file.openDataSet("Model/Elements/Number_of_Output_Fields");
	int Element_No_of_Output_Fields[Number_of_Elements[0]+1]; Element_No_of_Output_Fields_DataSet.read( Element_No_of_Output_Fields, PredType::NATIVE_INT);

	DataSet Element_Number_of_Gauss_Points_DataSet = file.openDataSet("Model/Elements/Number_of_Gauss_Points");
	int Element_Number_of_Gauss_Points[Number_of_Elements[0]+1]; Element_Number_of_Gauss_Points_DataSet.read( Element_Number_of_Gauss_Points, PredType::NATIVE_INT);

	hsize_t dims_out[2]; int ndims; DataSet Element_Outputs_DataSet = file.openDataSet("Model/Elements/Outputs");
	DataSpace Element_Outputs_DataSpace = Element_Outputs_DataSet.getSpace(); ndims = Element_Outputs_DataSpace.getSimpleExtentDims( dims_out, NULL);
	float Element_Outputs[dims_out[0]][dims_out[1]]; Element_Outputs_DataSet.read( Element_Outputs, PredType::NATIVE_FLOAT);

	///////////////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	Incremental_Energy_Database = new float [(Number_of_Elements[0]+1)*No_of_TimeSteps[0]];
  	Total_Energy_Database = new float [(Number_of_Elements[0]+1)*No_of_TimeSteps[0]];

	for(int Ctime=0; Ctime < No_of_TimeSteps[0];  Ctime++){

		for (int i = 0; i <= Number_of_Elements[0]; i++){

			int connectivity_index = Element_Index_to_Connectivity[i];

			if(connectivity_index!=-1){
		
				int NOF_Gauss_Points = Element_Number_of_Gauss_Points[i];
				int Output_Index = Element_Index_to_Outputs[i];
				int NOF_Output_Fields = Element_No_of_Output_Fields[i];
				float Energy_Increment=0;

				for(int j=0; j< NOF_Gauss_Points ; j++){
					float El_Strain_Tuple[9] ={Element_Outputs[Output_Index][Ctime],Element_Outputs[Output_Index+1][Ctime],Element_Outputs[Output_Index+3][Ctime],Element_Outputs[Output_Index+1][Ctime],Element_Outputs[Output_Index+2][Ctime],Element_Outputs[Output_Index+4][Ctime],Element_Outputs[Output_Index+3][Ctime],Element_Outputs[Output_Index+4][Ctime],Element_Outputs[Output_Index+5][Ctime]};
					Output_Index = Output_Index+6;
					float Pl_Strain_Tuple[9] ={Element_Outputs[Output_Index][Ctime],Element_Outputs[Output_Index+1][Ctime],Element_Outputs[Output_Index+3][Ctime],Element_Outputs[Output_Index+1][Ctime],Element_Outputs[Output_Index+2][Ctime],Element_Outputs[Output_Index+4][Ctime],Element_Outputs[Output_Index+3][Ctime],Element_Outputs[Output_Index+4][Ctime],Element_Outputs[Output_Index+5][Ctime]};
					Output_Index = Output_Index+6;
					float Stress_Tuple[9] ={Element_Outputs[Output_Index][Ctime],Element_Outputs[Output_Index+1][Ctime],Element_Outputs[Output_Index+3][Ctime],Element_Outputs[Output_Index+1][Ctime],Element_Outputs[Output_Index+2][Ctime],Element_Outputs[Output_Index+4][Ctime],Element_Outputs[Output_Index+3][Ctime],Element_Outputs[Output_Index+4][Ctime],Element_Outputs[Output_Index+5][Ctime]};
					Output_Index = Output_Index+6;

					Output_Index = Output_Index-18;

					if(Ctime>0){
						float Prev_El_Strain_Tuple[9] ={Element_Outputs[Output_Index][Ctime-1],Element_Outputs[Output_Index+1][Ctime-1],Element_Outputs[Output_Index+3][Ctime-1],Element_Outputs[Output_Index+1][Ctime-1],Element_Outputs[Output_Index+2][Ctime-1],Element_Outputs[Output_Index+4][Ctime-1],Element_Outputs[Output_Index+3][Ctime-1],Element_Outputs[Output_Index+4][Ctime-1],Element_Outputs[Output_Index+5][Ctime-1]};
						Output_Index = Output_Index+6;
						float Prev_Pl_Strain_Tuple[9] ={Element_Outputs[Output_Index][Ctime-1],Element_Outputs[Output_Index+1][Ctime-1],Element_Outputs[Output_Index+3][Ctime-1],Element_Outputs[Output_Index+1][Ctime-1],Element_Outputs[Output_Index+2][Ctime-1],Element_Outputs[Output_Index+4][Ctime-1],Element_Outputs[Output_Index+3][Ctime-1],Element_Outputs[Output_Index+4][Ctime-1],Element_Outputs[Output_Index+5][Ctime-1]};
						Output_Index = Output_Index+6;
						float Prev_Stress_Tuple[9] ={Element_Outputs[Output_Index][Ctime-1],Element_Outputs[Output_Index+1][Ctime-1],Element_Outputs[Output_Index+3][Ctime-1],Element_Outputs[Output_Index+1][Ctime-1],Element_Outputs[Output_Index+2][Ctime-1],Element_Outputs[Output_Index+4][Ctime-1],Element_Outputs[Output_Index+3][Ctime-1],Element_Outputs[Output_Index+4][Ctime-1],Element_Outputs[Output_Index+5][Ctime-1]};
						Output_Index = Output_Index+6;
						Energy_Increment+=(Stress_Tuple[0]-Prev_Stress_Tuple[0])*(El_Strain_Tuple[0]+Pl_Strain_Tuple[0]-Prev_El_Strain_Tuple[0]-Prev_Pl_Strain_Tuple[0])+
										(Stress_Tuple[1]-Prev_Stress_Tuple[1])*(El_Strain_Tuple[1]+Pl_Strain_Tuple[1]-Prev_El_Strain_Tuple[1]-Prev_Pl_Strain_Tuple[1])+
										(Stress_Tuple[2]-Prev_Stress_Tuple[2])*(El_Strain_Tuple[2]+Pl_Strain_Tuple[2]-Prev_El_Strain_Tuple[2]-Prev_Pl_Strain_Tuple[2])+
										(Stress_Tuple[3]-Prev_Stress_Tuple[3])*(El_Strain_Tuple[3]+Pl_Strain_Tuple[3]-Prev_El_Strain_Tuple[3]-Prev_Pl_Strain_Tuple[3])+
										(Stress_Tuple[4]-Prev_Stress_Tuple[4])*(El_Strain_Tuple[4]+Pl_Strain_Tuple[4]-Prev_El_Strain_Tuple[4]-Prev_Pl_Strain_Tuple[4])+
										(Stress_Tuple[5]-Prev_Stress_Tuple[5])*(El_Strain_Tuple[5]+Pl_Strain_Tuple[5]-Prev_El_Strain_Tuple[5]-Prev_Pl_Strain_Tuple[5])+
										(Stress_Tuple[6]-Prev_Stress_Tuple[6])*(El_Strain_Tuple[6]+Pl_Strain_Tuple[6]-Prev_El_Strain_Tuple[6]-Prev_Pl_Strain_Tuple[6])+
										(Stress_Tuple[7]-Prev_Stress_Tuple[7])*(El_Strain_Tuple[7]+Pl_Strain_Tuple[7]-Prev_El_Strain_Tuple[7]-Prev_Pl_Strain_Tuple[7])+
										(Stress_Tuple[8]-Prev_Stress_Tuple[8])*(El_Strain_Tuple[8]+Pl_Strain_Tuple[8]-Prev_El_Strain_Tuple[8]-Prev_Pl_Strain_Tuple[8]);
					}
					else{
						Energy_Increment+=(Stress_Tuple[0])*(El_Strain_Tuple[0]+Pl_Strain_Tuple[0])+
										(Stress_Tuple[1])*(El_Strain_Tuple[1]+Pl_Strain_Tuple[1])+
										(Stress_Tuple[2])*(El_Strain_Tuple[2]+Pl_Strain_Tuple[2])+
										(Stress_Tuple[3])*(El_Strain_Tuple[3]+Pl_Strain_Tuple[3])+
										(Stress_Tuple[4])*(El_Strain_Tuple[4]+Pl_Strain_Tuple[4])+
										(Stress_Tuple[5])*(El_Strain_Tuple[5]+Pl_Strain_Tuple[5])+
										(Stress_Tuple[6])*(El_Strain_Tuple[6]+Pl_Strain_Tuple[6])+
										(Stress_Tuple[7])*(El_Strain_Tuple[7]+Pl_Strain_Tuple[7])+
										(Stress_Tuple[8])*(El_Strain_Tuple[8]+Pl_Strain_Tuple[8]);
						Output_Index = Output_Index+18;
					}
				}

				Incremental_Energy_Database[Ctime*(Number_of_Elements[0]+1)+i] = Energy_Increment/8;
				if(Ctime>0)
					Total_Energy_Database[Ctime*(Number_of_Elements[0]+1)+i] = Total_Energy_Database[(Ctime-1)*(Number_of_Elements[0]+1)+i] + Incremental_Energy_Database[Ctime*(Number_of_Elements[0]+1)+i];
				else
					Total_Energy_Database[Ctime*(Number_of_Elements[0]+1)+i] = Incremental_Energy_Database[Ctime*(Number_of_Elements[0]+1)+i];
			}
		}
	}
	Energy_Database_Status=1;
}


void pvESSI::GetMesh(){

	/////////////////////////////////////////// Reading a HDF5 file /////////////////////////////////////////////////////////////////////////////////////////////////

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	//////////////////////////////////////// Reading General Outside Data ////////////////////////////////////////////////////////////////////////////////////////////

	DataSet No_of_Elements_DataSet = file.openDataSet("Number_of_Elements");
	int Number_of_Elements[1]; No_of_Elements_DataSet.read( Number_of_Elements, PredType::NATIVE_INT);

	DataSet Number_of_Nodes_DataSet = file.openDataSet("Number_of_Nodes");
	int Number_of_Nodes[1]; Number_of_Nodes_DataSet.read( Number_of_Nodes, PredType::NATIVE_INT);

	DataSet No_of_TimeSteps_DataSet = file.openDataSet("Number_of_Time_Steps");
	int No_of_TimeSteps[1]; No_of_TimeSteps_DataSet.read( No_of_TimeSteps, PredType::NATIVE_INT);

	//////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet Element_Index_to_Connectivity_DataSet = file.openDataSet("Model/Elements/Index_to_Connectivity");
	int *Element_Index_to_Connectivity; Element_Index_to_Connectivity = (int*) malloc((Number_of_Elements[0]+1) * sizeof(int)); Element_Index_to_Connectivity_DataSet.read( Element_Index_to_Connectivity, PredType::NATIVE_INT);

	DataSet Element_Index_to_Gauss_Point_Coordinates_DataSet = file.openDataSet("Model/Elements/Index_to_Gauss_Point_Coordinates");
	int *Element_Index_to_Gauss_Point_Coordinates; Element_Index_to_Gauss_Point_Coordinates = (int*) malloc((Number_of_Elements[0]+1) * sizeof(int)); Element_Index_to_Gauss_Point_Coordinates_DataSet.read( Element_Index_to_Gauss_Point_Coordinates, PredType::NATIVE_INT);

	DataSet Element_Number_of_Gauss_Points_DataSet = file.openDataSet("Model/Elements/Number_of_Gauss_Points");
	int *Element_Number_of_Gauss_Points; Element_Number_of_Gauss_Points = (int*) malloc((Number_of_Elements[0]+1) * sizeof(int)); Element_Number_of_Gauss_Points_DataSet.read( Element_Number_of_Gauss_Points, PredType::NATIVE_INT);

	DataSet Element_Number_of_Nodes_DataSet = file.openDataSet("Model/Elements/Number_of_Nodes");
	int *Element_Number_of_Nodes; Element_Number_of_Nodes = (int*) malloc((Number_of_Elements[0]+1) * sizeof(int)); Element_Number_of_Nodes_DataSet.read( Element_Number_of_Nodes, PredType::NATIVE_INT);

	DataSet Element_Connectivity_DataSet = file.openDataSet("Model/Elements/Connectivity");
	DataSpace Element_Connectivity_DataSpace = Element_Connectivity_DataSet.getSpace(); hsize_t dims_out[1]; int ndims = Element_Connectivity_DataSpace.getSimpleExtentDims( dims_out, NULL);
	int *Element_Connectivity; Element_Connectivity = (int*) malloc(dims_out[0] * sizeof(int));	Element_Connectivity_DataSet.read( Element_Connectivity, PredType::NATIVE_INT);

	DataSet Element_Gauss_Point_Coordinates_DataSet = file.openDataSet("Model/Elements/Gauss_Point_Coordinates");
	DataSpace Element_Gauss_Point_Coordinates_DataSpace = Element_Gauss_Point_Coordinates_DataSet.getSpace(); ndims = Element_Gauss_Point_Coordinates_DataSpace.getSimpleExtentDims( dims_out, NULL);
	float *Element_Gauss_Point_Coordinates; Element_Gauss_Point_Coordinates = (float*) malloc(dims_out[0] * sizeof(float));	Element_Gauss_Point_Coordinates_DataSet.read( Element_Gauss_Point_Coordinates, PredType::NATIVE_INT);

	///////////////////////////////////////////////// Reading Nodes datat //////////////////////////////////////////////////////////////////////////////////////////////

	DataSet Node_Coordinates_DataSet = file.openDataSet("Model/Nodes/Coordinates");
	float *Node_Coordinates; Node_Coordinates = (float*) malloc((Number_of_Nodes[0]*3) * sizeof(float)); Node_Coordinates_DataSet.read( Node_Coordinates, PredType::NATIVE_FLOAT);

	DataSet Node_Index_to_Coordinates_DataSet = file.openDataSet("Model/Nodes/Index_to_Coordinates");
	int *Node_Index_to_Coordinates; Node_Index_to_Coordinates = (int*) malloc((Number_of_Nodes[0]+1) * sizeof(int)); Node_Index_to_Coordinates_DataSet.read( Node_Index_to_Coordinates, PredType::NATIVE_INT);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	points->SetNumberOfPoints(Number_of_Nodes[0]);
	UGrid_Mesh->Allocate(Number_of_Elements[0]);

	for (int i = 0; i <= Number_of_Nodes[0]; i++){
		if (Node_Index_to_Coordinates[i]!=-1){
			points->InsertPoint(i, Node_Coordinates[Node_Index_to_Coordinates[i]],Node_Coordinates[Node_Index_to_Coordinates[i]+1],Node_Coordinates[Node_Index_to_Coordinates[i]+2]);
		}
	}

	UGrid_Mesh->SetPoints(points);

	/////////////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	for (int i = 0; i <= Number_of_Elements[0]; i++){

		int connectivity_index = Element_Index_to_Connectivity[i];

		if(connectivity_index!=-1){
	
			int No_of_Element_Nodes = Element_Number_of_Nodes[i];
			vtkIdType Vertices[No_of_Element_Nodes];

			int Cell_Type = ESSI_to_VTK_Element.find(No_of_Element_Nodes)->second;
			std::vector<int> Nodes_Connectivity_Order = ESSI_to_VTK_Connectivity.find(No_of_Element_Nodes)->second;

			for(int j=0; j<No_of_Element_Nodes ; j++){
				Vertices[j] = Element_Connectivity[connectivity_index+Nodes_Connectivity_Order[j]];
			}

			UGrid_Mesh->InsertNextCell(Cell_Type, No_of_Element_Nodes, Vertices);

			// Material_Properties->InsertValue(i,i);
		}
	}

	mesh_build=1;
	
}

void pvESSI::SetMetaArrays( vtkSmartPointer<vtkFloatArray> &vtk_Generalized_Displacements, vtkSmartPointer<vtkFloatArray> &vtk_Generalized_Velocity, vtkSmartPointer<vtkFloatArray> &vtk_Generalized_Acceleration, vtkSmartPointer<vtkFloatArray> &Elastic_Strain_Tensor, vtkSmartPointer<vtkFloatArray> &Plastic_Strain_Tensor, vtkSmartPointer<vtkFloatArray> &Stress_Tensor, vtkSmartPointer<vtkFloatArray> &Material_Properties, vtkSmartPointer<vtkFloatArray> &Total_Energy, vtkSmartPointer<vtkFloatArray> &Incremental_Energy){

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

	/*************************************************************************************************************************************************/
}