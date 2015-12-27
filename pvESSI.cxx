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
	this->SetNumberOfOutputPorts(2);
	Energy_Database_Status=-1;
	UGrid_Mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
	UGrid_Gauss_Mesh = vtkSmartPointer<vtkUnstructuredGrid>:: New();
	UGrid_All_Mesh = new vtkSmartPointer<vtkUnstructuredGrid>[2];
	this->set_VTK_To_ESSI_Elements_Connectivity();
	Current_Time =0;
}

int pvESSI::RequestData(vtkInformation *vtkNotUsed(request),vtkInformationVector **vtkNotUsed(inputVector),	vtkInformationVector *outputVector){
 
	// get the info object
	vtkInformation *outInfo1 = outputVector->GetInformationObject(0);
	vtkInformation *outInfo2 = outputVector->GetInformationObject(1);
	// outInfo->Print(std::cout);

	if (!whetehr_general_mesh_build){ this->GetGaussMesh(); this->GetGeneralMesh();} 
	
	// if(Energy_Database_Status==-1)
	// 	this->Make_Energy_Database();

	// int extent[6] = {0,-1,0,-1,0,-1};
	// outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);

  	Current_Time = outInfo1->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
  	// int Ctime2 = outInfo2->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
  	// int Clength = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  	// double* Csteps = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

 	cout<< "Current_Time " << Current_Time<< endl;
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
	vtkUnstructuredGrid *GeneralMesh = vtkUnstructuredGrid::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));
	vtkUnstructuredGrid *GaussMesh = vtkUnstructuredGrid::SafeDownCast(outInfo2->Get(vtkDataObject::DATA_OBJECT()));

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


	Build_Gauss_Attributes();
	Build_Node_Attributes();


	GeneralMesh->ShallowCopy(UGrid_Mesh);
	GaussMesh->ShallowCopy(UGrid_Gauss_Mesh);

	// output->ShallowCopy(UGrid_All_Mesh);

	return 1;
}

int pvESSI::RequestInformation( vtkInformation *request, vtkInformationVector **vtkNotUsed(inVec), vtkInformationVector* outVec){

	vtkInformation* outInfo1 = outVec->GetInformationObject(0);
	vtkInformation* outInfo2 = outVec->GetInformationObject(1);

	//////////////////////////////////////////////////////// Reading Hdf5 File ///////////////////////////////////////////////////////////////////////////////////////

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	////////////////////////////////////////////////////// Reading General Outside Data /////////////////////////////////////////////////////////////////////////////

	DataSet No_of_Elements_DataSet = file.openDataSet("Number_of_Elements");
	int Number_of_Elements[1]; No_of_Elements_DataSet.read( Number_of_Elements, PredType::NATIVE_INT);

	DataSet Number_of_Nodes_DataSet = file.openDataSet("Number_of_Nodes");
	int Number_of_Nodes[1]; Number_of_Nodes_DataSet.read( Number_of_Nodes, PredType::NATIVE_INT);

	DataSet No_of_TimeSteps_DataSet = file.openDataSet("Number_of_Time_Steps");
	int No_of_TimeSteps[1]; No_of_TimeSteps_DataSet.read( No_of_TimeSteps, PredType::NATIVE_INT);

	DataSet Element_Number_of_Gauss_Points_DataSet = file.openDataSet("Model/Elements/Number_of_Gauss_Points");
	int Element_Number_of_Gauss_Points[Number_of_Elements[0]+1]; Element_Number_of_Gauss_Points_DataSet.read( Element_Number_of_Gauss_Points, PredType::NATIVE_INT);

	this->Number_of_Elements = Number_of_Elements[0];
	this->Number_of_Nodes = Number_of_Nodes[0];
	this->Number_of_Gauss_Nodes = 0;

	///////////////////////////////////////////////////////////////// Evaluating Number of Gauss Points ////////////////////////////////////////////////////////////////////////////////////////////

    for (int i=0; i <= this->Number_of_Elements; i++){
        if (Element_Number_of_Gauss_Points[i] > 0){
            this->Number_of_Gauss_Nodes += Element_Number_of_Gauss_Points[i];
        }
    }

    float *Prev_Total_Energy_Database; Prev_Total_Energy_Database = (float*) malloc((Number_of_Gauss_Nodes)* sizeof(float)); 
    for (int i=0; i < this->Number_of_Gauss_Nodes; i++){
    	Prev_Total_Energy_Database[i]=0;
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	this->Number_Of_Time_Steps = No_of_TimeSteps[0];
	double Time_Steps[No_of_TimeSteps[0]]; Time_Steps[0]=0;

	for (int i =0;i<No_of_TimeSteps[0] ;i++)
		Time_Steps[i] = i;

	double Time_range[2]={0,No_of_TimeSteps[0]};

	// int extent[6]={0,100,0,100,0,100}; 


	// outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),extent,6);
	outInfo1->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),Time_Steps, No_of_TimeSteps[0]);
	outInfo1->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),Time_range,2);

	outInfo2->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),Time_Steps, No_of_TimeSteps[0]);
	outInfo2->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),Time_range,2);
	// outInfo->Set(vtkAlgorithm::CAN_PRODUCE_SUB_EXTENT(),1);


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

	ESSI_to_VTK_Element[1] = VTK_VERTEX;		ESSI_to_VTK_Element[2]  = VTK_LINE;					ESSI_to_VTK_Element[4] = VTK_QUAD;
	ESSI_to_VTK_Element[8] = VTK_HEXAHEDRON;	ESSI_to_VTK_Element[20] = VTK_QUADRATIC_HEXAHEDRON;	ESSI_to_VTK_Element[27] = VTK_TRIQUADRATIC_HEXAHEDRON;

	int Node[1] = {0};						connectivity_vector.assign(Node,Node+1);				ESSI_to_VTK_Connectivity[1] = connectivity_vector;	connectivity_vector.clear();
	int Line[2] = {0,1}; 					connectivity_vector.assign(Line,Line+2); 				ESSI_to_VTK_Connectivity[2] = connectivity_vector;	connectivity_vector.clear();
	int Quadrangle[4] = {0,1,2,3}; 			connectivity_vector.assign(Quadrangle,Quadrangle+4); 	ESSI_to_VTK_Connectivity[4] = connectivity_vector;	connectivity_vector.clear();
	int Hexahedron[8] = {4,5,6,7,0,1,2,3}; 	connectivity_vector.assign(Hexahedron,Hexahedron+8);	ESSI_to_VTK_Connectivity[8] = connectivity_vector;	connectivity_vector.clear();
	int Qudratic_hexahedron[20] = {6,5,4,7,2,1,0,3,13,12,15,14,9,8,11,10,18,17,16,19}; 	connectivity_vector.assign(Qudratic_hexahedron,Qudratic_hexahedron+20);	ESSI_to_VTK_Connectivity[20]= connectivity_vector;	connectivity_vector.clear();
	int TriQudratic_hexahedron[27] = {6,5,4,7,2,1,0,3,13,12,15,14,9,8,11,10,18,17,16,19,23,21,22,24,26,25,20}; 	connectivity_vector.assign(TriQudratic_hexahedron,TriQudratic_hexahedron+27);	ESSI_to_VTK_Connectivity[27]= connectivity_vector;	connectivity_vector.clear();
}

void pvESSI::Build_Node_Attributes(){

	/////////////////////////////////////////////////////////////////// Reading a HDF5 file /////////////////////////////////////////////////////////////////////////////////////////////////

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	DataSet Node_Index_to_Coordinates_DataSet = file.openDataSet("Model/Nodes/Index_to_Coordinates");
	int *Node_Index_to_Coordinates; Node_Index_to_Coordinates = (int*) malloc((Number_of_Nodes+1) * sizeof(int)); Node_Index_to_Coordinates_DataSet.read( Node_Index_to_Coordinates, PredType::NATIVE_INT);

	DataSet Node_Index_to_Generalized_Displacements_DataSet = file.openDataSet("Model/Nodes/Index_to_Generalized_Displacements");
	int *Node_Index_to_Generalized_Displacements; Node_Index_to_Generalized_Displacements = (int*) malloc((Number_of_Nodes+1) * sizeof(int)); Node_Index_to_Generalized_Displacements_DataSet.read( Node_Index_to_Generalized_Displacements, PredType::NATIVE_INT);

	///////////////////////////////////////////  Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	
	
    DataSet Node_Generalized_Displacements_DataSet = file.openDataSet("Model/Nodes/Generalized_Displacements"); hsize_t dims_out[2]; int ndims;
	DataSpace Node_Generalized_Displacements_DataSpace = Node_Generalized_Displacements_DataSet.getSpace(); ndims = Node_Generalized_Displacements_DataSpace.getSimpleExtentDims( dims_out, NULL);
	hsize_t col_dims[1];  col_dims[0] = dims_out[0]; DataSpace mspace2( 1, col_dims ); 
	hsize_t offset[2] = { 0, Current_Time }; hsize_t  count[2] = { dims_out[0], 1 }; /* Define the column (hyperslab) to read. */
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

	for (int i = 0; i <= Number_of_Nodes; i++){

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

	/////////////////////////////////////////////////////////////////////// Enable the Data Array  to show in Paraview  ///////////////////////////////////////////////////////////////

	UGrid_Mesh->GetPointData()->AddArray(vtk_Generalized_Displacements);
  	// UGrid_Mesh->GetPointData()->AddArray(vtk_Generalized_Velocity);
  	// UGrid_Mesh->GetPointData()->AddArray(vtk_Generalized_Acceleration);
}

void pvESSI::Build_Gauss_Attributes(){

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	//////////////////////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet Element_Index_to_Outputs_DataSet = file.openDataSet("Model/Elements/Index_to_Outputs");
	int Element_Index_to_Outputs[Number_of_Elements+1]; Element_Index_to_Outputs_DataSet.read( Element_Index_to_Outputs, PredType::NATIVE_INT);

	DataSet Element_Number_of_Gauss_Points_DataSet = file.openDataSet("Model/Elements/Number_of_Gauss_Points");
	int Element_Number_of_Gauss_Points[Number_of_Elements+1]; Element_Number_of_Gauss_Points_DataSet.read( Element_Number_of_Gauss_Points, PredType::NATIVE_INT);

	DataSet Element_Outputs_DataSet = file.openDataSet("Model/Elements/Outputs"); hsize_t dims_out[2]; int ndims;
	DataSpace Element_Outputs_DataSpace = Element_Outputs_DataSet.getSpace(); ndims = Element_Outputs_DataSpace.getSimpleExtentDims( dims_out, NULL);
	hsize_t col_dims[1];  col_dims[0] = dims_out[0]; DataSpace mspace1( 1, col_dims );
    hsize_t offset[2] = { 0, Current_Time }; hsize_t  count[2] = { dims_out[0], 1 }; /* Define the column (hyperslab) to read. */
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

 	SetMetaArrays( vtk_Generalized_Displacements, vtk_Generalized_Velocity, vtk_Generalized_Acceleration, Elastic_Strain_Tensor, Plastic_Strain_Tensor, Stress_Tensor, Material_Properties, Total_Energy, Incremental_Energy);
 	int Index_of_Gauss_Nodes =0;

	for (int i = 0; i <= Number_of_Elements; i++){

		if (Element_Number_of_Gauss_Points[i] > 0){

			/************************************************** Gauss Calculations **********************************************/

			int No_of_Element_Gauss_Nodes = Element_Number_of_Gauss_Points[i];
			int Output_Index = Element_Index_to_Outputs[i];

			for(int j=0; j< No_of_Element_Gauss_Nodes ; j++){

				float Gauss_Incremental_Energy[1] = {0}, Gauss_Total_Energy[1]= {0} ;
				float El_Strain_Tuple[9] ={Element_Outputs[Output_Index],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+2],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+5]};
				Output_Index = Output_Index+6;
				float Pl_Strain_Tuple[9] ={Element_Outputs[Output_Index],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+2],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+5]};
				Output_Index = Output_Index+6;
				float Stress_Tuple[9] ={Element_Outputs[Output_Index],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+2],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+5]};
				Output_Index = Output_Index+6;
				Elastic_Strain_Tensor->InsertTupleValue(Index_of_Gauss_Nodes,El_Strain_Tuple);
				Plastic_Strain_Tensor->InsertTupleValue(Index_of_Gauss_Nodes,Pl_Strain_Tuple);
				Stress_Tensor->InsertTupleValue(Index_of_Gauss_Nodes,Stress_Tuple);

				if(Current_Time>0){

					DataSet Element_Outputs_DataSet = file.openDataSet("Model/Elements/Outputs"); hsize_t dims_out[2]; int ndims;
					DataSpace Element_Outputs_DataSpace = Element_Outputs_DataSet.getSpace(); ndims = Element_Outputs_DataSpace.getSimpleExtentDims( dims_out, NULL);
					hsize_t col_dims[1];  col_dims[0] = dims_out[0]; DataSpace mspace1( 1, col_dims );
				    hsize_t offset[2] = { 0, Current_Time-1 }; hsize_t  count[2] = { dims_out[0], 1 }; /* Define the column (hyperslab) to read. */
				    float *Element_Outputs; Element_Outputs = (float*) malloc(dims_out[0] * sizeof(float)); // buffer for column to be read and defining hyperslab and read.
				    Element_Outputs_DataSpace.selectHyperslab( H5S_SELECT_SET, count, offset ); Element_Outputs_DataSet.read( Element_Outputs, PredType::NATIVE_FLOAT, mspace1, Element_Outputs_DataSpace );

					float Prev_El_Strain_Tuple[9] ={Element_Outputs[Output_Index],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+2],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+5]};
					Output_Index = Output_Index+6;
					float Prev_Pl_Strain_Tuple[9] ={Element_Outputs[Output_Index],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+2],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+5]};
					Output_Index = Output_Index+6;
					float Prev_Stress_Tuple[9] ={Element_Outputs[Output_Index],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+1],Element_Outputs[Output_Index+2],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+3],Element_Outputs[Output_Index+4],Element_Outputs[Output_Index+5]};
					Output_Index = Output_Index+6;

					Gauss_Incremental_Energy[0]+=(Stress_Tuple[0]-Prev_Stress_Tuple[0])*(El_Strain_Tuple[0]+Pl_Strain_Tuple[0]-Prev_El_Strain_Tuple[0]-Prev_Pl_Strain_Tuple[0])+
											  (Stress_Tuple[1]-Prev_Stress_Tuple[1])*(El_Strain_Tuple[1]+Pl_Strain_Tuple[1]-Prev_El_Strain_Tuple[1]-Prev_Pl_Strain_Tuple[1])+
											  (Stress_Tuple[2]-Prev_Stress_Tuple[2])*(El_Strain_Tuple[2]+Pl_Strain_Tuple[2]-Prev_El_Strain_Tuple[2]-Prev_Pl_Strain_Tuple[2])+
											  (Stress_Tuple[3]-Prev_Stress_Tuple[3])*(El_Strain_Tuple[3]+Pl_Strain_Tuple[3]-Prev_El_Strain_Tuple[3]-Prev_Pl_Strain_Tuple[3])+
											  (Stress_Tuple[4]-Prev_Stress_Tuple[4])*(El_Strain_Tuple[4]+Pl_Strain_Tuple[4]-Prev_El_Strain_Tuple[4]-Prev_Pl_Strain_Tuple[4])+
											  (Stress_Tuple[5]-Prev_Stress_Tuple[5])*(El_Strain_Tuple[5]+Pl_Strain_Tuple[5]-Prev_El_Strain_Tuple[5]-Prev_Pl_Strain_Tuple[5])+
											  (Stress_Tuple[6]-Prev_Stress_Tuple[6])*(El_Strain_Tuple[6]+Pl_Strain_Tuple[6]-Prev_El_Strain_Tuple[6]-Prev_Pl_Strain_Tuple[6])+
											  (Stress_Tuple[7]-Prev_Stress_Tuple[7])*(El_Strain_Tuple[7]+Pl_Strain_Tuple[7]-Prev_El_Strain_Tuple[7]-Prev_Pl_Strain_Tuple[7])+
											  (Stress_Tuple[8]-Prev_Stress_Tuple[8])*(El_Strain_Tuple[8]+Pl_Strain_Tuple[8]-Prev_El_Strain_Tuple[8]-Prev_Pl_Strain_Tuple[8]);
					
					Gauss_Total_Energy[0] = Prev_Total_Energy_Database[Index_of_Gauss_Nodes] + Gauss_Incremental_Energy[0];
			
				}
				else{

					Gauss_Incremental_Energy[0]=(Stress_Tuple[0])*(El_Strain_Tuple[0]+Pl_Strain_Tuple[0])+
			  			 					  	(Stress_Tuple[1])*(El_Strain_Tuple[1]+Pl_Strain_Tuple[1])+
			  			 					  	(Stress_Tuple[2])*(El_Strain_Tuple[2]+Pl_Strain_Tuple[2])+
			  			 					  	(Stress_Tuple[3])*(El_Strain_Tuple[3]+Pl_Strain_Tuple[3])+
			  			 					  	(Stress_Tuple[4])*(El_Strain_Tuple[4]+Pl_Strain_Tuple[4])+
			  			 					  	(Stress_Tuple[5])*(El_Strain_Tuple[5]+Pl_Strain_Tuple[5])+
			  			 					  	(Stress_Tuple[6])*(El_Strain_Tuple[6]+Pl_Strain_Tuple[6])+
			  			 					  	(Stress_Tuple[7])*(El_Strain_Tuple[7]+Pl_Strain_Tuple[7])+
			  			 					  	(Stress_Tuple[8])*(El_Strain_Tuple[8]+Pl_Strain_Tuple[8]);

					Gauss_Total_Energy[0] = Gauss_Incremental_Energy[0];
					Prev_Total_Energy_Database[Index_of_Gauss_Nodes] = Gauss_Total_Energy[0];
				}

				Incremental_Energy->InsertTupleValue(Index_of_Gauss_Nodes,Gauss_Incremental_Energy);
				Total_Energy->InsertTupleValue(Index_of_Gauss_Nodes,Gauss_Total_Energy);

			}

			Index_of_Gauss_Nodes+=1;

			/********************************************************************************************************************/
		}
	}

	UGrid_Gauss_Mesh->GetCellData()->AddArray(Stress_Tensor);
  	UGrid_Gauss_Mesh->GetCellData()->AddArray(Plastic_Strain_Tensor);
  	UGrid_Gauss_Mesh->GetCellData()->AddArray(Elastic_Strain_Tensor);
  	UGrid_Gauss_Mesh->GetCellData()->AddArray(Total_Energy);
	UGrid_Gauss_Mesh->GetCellData()->AddArray(Incremental_Energy);

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
}


void pvESSI::GetGeneralMesh(){

	/////////////////////////////////////////// Reading a HDF5 file /////////////////////////////////////////////////////////////////////////////////////////////////

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	//////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet Element_Index_to_Connectivity_DataSet = file.openDataSet("Model/Elements/Index_to_Connectivity");
	int *Element_Index_to_Connectivity; Element_Index_to_Connectivity = (int*) malloc((Number_of_Elements+1) * sizeof(int)); Element_Index_to_Connectivity_DataSet.read( Element_Index_to_Connectivity, PredType::NATIVE_INT);

	DataSet Element_Index_to_Gauss_Point_Coordinates_DataSet = file.openDataSet("Model/Elements/Index_to_Gauss_Point_Coordinates");
	int *Element_Index_to_Gauss_Point_Coordinates; Element_Index_to_Gauss_Point_Coordinates = (int*) malloc((Number_of_Elements+1) * sizeof(int)); Element_Index_to_Gauss_Point_Coordinates_DataSet.read( Element_Index_to_Gauss_Point_Coordinates, PredType::NATIVE_INT);

	DataSet Element_Number_of_Gauss_Points_DataSet = file.openDataSet("Model/Elements/Number_of_Gauss_Points");
	int *Element_Number_of_Gauss_Points; Element_Number_of_Gauss_Points = (int*) malloc((Number_of_Elements+1) * sizeof(int)); Element_Number_of_Gauss_Points_DataSet.read( Element_Number_of_Gauss_Points, PredType::NATIVE_INT);

	DataSet Element_Number_of_Nodes_DataSet = file.openDataSet("Model/Elements/Number_of_Nodes");
	int *Element_Number_of_Nodes; Element_Number_of_Nodes = (int*) malloc((Number_of_Elements+1) * sizeof(int)); Element_Number_of_Nodes_DataSet.read( Element_Number_of_Nodes, PredType::NATIVE_INT);

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
	int *Node_Index_to_Coordinates; Node_Index_to_Coordinates = (int*) malloc((Number_of_Nodes+1) * sizeof(int)); Node_Index_to_Coordinates_DataSet.read( Node_Index_to_Coordinates, PredType::NATIVE_INT);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	points->SetNumberOfPoints(Number_of_Nodes);
	UGrid_Mesh->Allocate(Number_of_Elements);

	for (int i = 0; i <= Number_of_Nodes; i++){
		if (Node_Index_to_Coordinates[i]!=-1){
			points->InsertPoint(i, Node_Coordinates[Node_Index_to_Coordinates[i]],Node_Coordinates[Node_Index_to_Coordinates[i]+1],Node_Coordinates[Node_Index_to_Coordinates[i]+2]);
		}
	}

	UGrid_Mesh->SetPoints(points);

	/////////////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	for (int i = 0; i <= Number_of_Elements; i++){

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

		}
	}

	whetehr_general_mesh_build=1;
	
}

void pvESSI::GetGaussMesh(){

	/////////////////////////////////////////// Reading a HDF5 file /////////////////////////////////////////////////////////////////////////////////////////////////

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	/////////////////////////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet Element_Index_to_Connectivity_DataSet = file.openDataSet("Model/Elements/Index_to_Connectivity");
	int *Element_Index_to_Connectivity; Element_Index_to_Connectivity = (int*) malloc((Number_of_Elements+1) * sizeof(int)); Element_Index_to_Connectivity_DataSet.read( Element_Index_to_Connectivity, PredType::NATIVE_INT);

	DataSet Element_Index_to_Gauss_Point_Coordinates_DataSet = file.openDataSet("Model/Elements/Index_to_Gauss_Point_Coordinates");
	int *Element_Index_to_Gauss_Point_Coordinates; Element_Index_to_Gauss_Point_Coordinates = (int*) malloc((Number_of_Elements+1) * sizeof(int)); Element_Index_to_Gauss_Point_Coordinates_DataSet.read( Element_Index_to_Gauss_Point_Coordinates, PredType::NATIVE_INT);

	DataSet Element_Number_of_Gauss_Points_DataSet = file.openDataSet("Model/Elements/Number_of_Gauss_Points");
	int *Element_Number_of_Gauss_Points; Element_Number_of_Gauss_Points = (int*) malloc((Number_of_Elements+1) * sizeof(int)); Element_Number_of_Gauss_Points_DataSet.read( Element_Number_of_Gauss_Points, PredType::NATIVE_INT);

	DataSet Element_Connectivity_DataSet = file.openDataSet("Model/Elements/Connectivity");
	DataSpace Element_Connectivity_DataSpace = Element_Connectivity_DataSet.getSpace(); hsize_t dims_out[1]; int ndims = Element_Connectivity_DataSpace.getSimpleExtentDims( dims_out, NULL);
	int *Element_Connectivity; Element_Connectivity = (int*) malloc(dims_out[0] * sizeof(int));	Element_Connectivity_DataSet.read( Element_Connectivity, PredType::NATIVE_INT);

	DataSet Element_Gauss_Point_Coordinates_DataSet = file.openDataSet("Model/Elements/Gauss_Point_Coordinates");
	DataSpace Element_Gauss_Point_Coordinates_DataSpace = Element_Gauss_Point_Coordinates_DataSet.getSpace(); ndims = Element_Gauss_Point_Coordinates_DataSpace.getSimpleExtentDims( dims_out, NULL);
	float *Element_Gauss_Point_Coordinates; Element_Gauss_Point_Coordinates = (float*) malloc(dims_out[0] * sizeof(float));	Element_Gauss_Point_Coordinates_DataSet.read( Element_Gauss_Point_Coordinates, PredType::NATIVE_FLOAT);

	////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(Number_of_Gauss_Nodes); int Index_of_Gauss_Nodes=0;

	for (int i=0; i <= Number_of_Elements; i++){
		if (Element_Number_of_Gauss_Points[i] > 0){
			int No_of_Element_Gauss_Nodes = Element_Number_of_Gauss_Points[i];
			for(int j=0; j<No_of_Element_Gauss_Nodes ; j++){
				points->InsertPoint(Index_of_Gauss_Nodes, Element_Gauss_Point_Coordinates[Element_Index_to_Gauss_Point_Coordinates[i]+3*j],Element_Gauss_Point_Coordinates[Element_Index_to_Gauss_Point_Coordinates[i]+3*j+1],Element_Gauss_Point_Coordinates[Element_Index_to_Gauss_Point_Coordinates[i]+3*j+2]);
				Index_of_Gauss_Nodes = Index_of_Gauss_Nodes+1;
			}
		}
	}

	UGrid_Gauss_Mesh->SetPoints(points);
	UGrid_Gauss_Mesh->Allocate(Number_of_Gauss_Nodes); vtkIdType OneVertex; Index_of_Gauss_Nodes=0;

	for (int i=0; i <= Number_of_Elements; i++){
		if (Element_Number_of_Gauss_Points[i] > 0){
			int No_of_Element_Gauss_Nodes = Element_Number_of_Gauss_Points[i];
			for(int j=0; j<No_of_Element_Gauss_Nodes ; j++){
				OneVertex =Index_of_Gauss_Nodes;
				UGrid_Gauss_Mesh->InsertNextCell(VTK_VERTEX, 1, &OneVertex);
				Index_of_Gauss_Nodes = Index_of_Gauss_Nodes+1;
			}
		}
	}

	whetehr_gauss_mesh_build=1;
	
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