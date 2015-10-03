#include "pvESSI.h"
 
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
	this->set_VTK_To_ESSI_Elements_Connectivity();	
	}
 
int pvESSI::RequestData(vtkInformation *vtkNotUsed(request),vtkInformationVector **vtkNotUsed(inputVector),	vtkInformationVector *outputVector){
 
	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	// outInfo->Print(std::cout);

	if(Energy_Database_Status==-1)
		this->Make_Energy_Database();

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

	// get the ouptut
	vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

		///////////////////////////////////////////// Reading a HDF5 file ////////////////////////////////////////////////////////////////////////////////////////////////////

	H5File file = H5File(this->FileName, H5F_ACC_RDONLY );

	//////////////////////////////////////// Reading General Outside Data ////////////////////////////////////////////////////////////////////////////////////////////

	DataSet No_of_Elements_DataSet = file.openDataSet("Number_of_Elements");
	int Number_of_Elements[1]; No_of_Elements_DataSet.read( Number_of_Elements, PredType::NATIVE_INT);

	DataSet Number_of_Nodes_DataSet = file.openDataSet("Number_of_Nodes");
	int Number_of_Nodes[1]; Number_of_Nodes_DataSet.read( Number_of_Nodes, PredType::NATIVE_INT);

	DataSet No_of_TimeSteps_DataSet = file.openDataSet("Number_of_Time_Steps");
	int No_of_TimeSteps[1]; No_of_TimeSteps_DataSet.read( No_of_TimeSteps, PredType::NATIVE_INT);

	DataSet Time_DataSet = file.openDataSet("time");
	int Time[No_of_TimeSteps[0]]; Time_DataSet.read( Time, PredType::NATIVE_INT);

	//////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet Element_Index_to_Connectivity_DataSet = file.openDataSet("Model/Elements/Index_to_Connectivity");
	int Element_Index_to_Connectivity[Number_of_Elements[0]+1]; Element_Index_to_Connectivity_DataSet.read( Element_Index_to_Connectivity, PredType::NATIVE_INT);

	DataSet Element_Index_to_Gauss_Point_Coordinates_DataSet = file.openDataSet("Model/Elements/Index_to_Gauss_Point_Coordinates");
	int Element_Index_to_Gauss_Point_Coordinates[Number_of_Elements[0]+1]; Element_Index_to_Gauss_Point_Coordinates_DataSet.read( Element_Index_to_Gauss_Point_Coordinates, PredType::NATIVE_INT);

	DataSet Element_Index_to_Outputs_DataSet = file.openDataSet("Model/Elements/Index_to_Outputs");
	int Element_Index_to_Outputs[Number_of_Elements[0]+1]; Element_Index_to_Outputs_DataSet.read( Element_Index_to_Outputs, PredType::NATIVE_INT);

	DataSet Element_No_of_Output_Fields_DataSet = file.openDataSet("Model/Elements/Number_of_Output_Fields");
	int Element_No_of_Output_Fields[Number_of_Elements[0]+1]; Element_No_of_Output_Fields_DataSet.read( Element_No_of_Output_Fields, PredType::NATIVE_INT);

	DataSet Element_Number_of_Gauss_Points_DataSet = file.openDataSet("Model/Elements/Number_of_Gauss_Points");
	int Element_Number_of_Gauss_Points[Number_of_Elements[0]+1]; Element_Number_of_Gauss_Points_DataSet.read( Element_Number_of_Gauss_Points, PredType::NATIVE_INT);

	DataSet Element_Number_of_Nodes_DataSet = file.openDataSet("Model/Elements/Number_of_Nodes");
	int Element_Number_of_Nodes[Number_of_Elements[0]+1]; Element_Number_of_Nodes_DataSet.read( Element_Number_of_Nodes, PredType::NATIVE_INT);

	DataSet Element_Connectivity_DataSet = file.openDataSet("Model/Elements/Connectivity");
	DataSpace Element_Connectivity_DataSpace = Element_Connectivity_DataSet.getSpace(); hsize_t dims_out[2]; int ndims = Element_Connectivity_DataSpace.getSimpleExtentDims( dims_out, NULL);
	int Element_Connectivity[dims_out[0]]; Element_Connectivity_DataSet.read( Element_Connectivity, PredType::NATIVE_INT);

	DataSet Element_Gauss_Point_Coordinates_DataSet = file.openDataSet("Model/Elements/Gauss_Point_Coordinates");
	DataSpace Element_Gauss_Point_Coordinates_DataSpace = Element_Gauss_Point_Coordinates_DataSet.getSpace(); ndims = Element_Gauss_Point_Coordinates_DataSpace.getSimpleExtentDims( dims_out, NULL);
	float Element_Gauss_Point_Coordinates[dims_out[0]]; Element_Gauss_Point_Coordinates_DataSet.read( Element_Gauss_Point_Coordinates, PredType::NATIVE_FLOAT);

	DataSet Element_Outputs_DataSet = file.openDataSet("Model/Elements/Outputs");
	DataSpace Element_Outputs_DataSpace = Element_Outputs_DataSet.getSpace(); ndims = Element_Outputs_DataSpace.getSimpleExtentDims( dims_out, NULL);
	float Element_Outputs[dims_out[0]][dims_out[1]]; Element_Outputs_DataSet.read( Element_Outputs, PredType::NATIVE_FLOAT);

	/////////////////////////////////////////////////// Reading Nodes datat //////////////////////////////////////////////////////////////////////////////////////////////

	DataSet Node_Coordinates_DataSet = file.openDataSet("Model/Nodes/Coordinates");
	float Node_Coordinates[Number_of_Nodes[0]*3]; Node_Coordinates_DataSet.read( Node_Coordinates, PredType::NATIVE_FLOAT);

	DataSet Node_Index_to_Coordinates_DataSet = file.openDataSet("Model/Nodes/Index_to_Coordinates");
	int Node_Index_to_Coordinates[Number_of_Nodes[0]+1]; Node_Index_to_Coordinates_DataSet.read( Node_Index_to_Coordinates, PredType::NATIVE_INT);

	DataSet Node_Generalized_Displacements_DataSet = file.openDataSet("Model/Nodes/Generalized_Displacements");
	DataSpace Node_Generalized_Displacements_DataSpace = Node_Generalized_Displacements_DataSet.getSpace(); ndims = Node_Generalized_Displacements_DataSpace.getSimpleExtentDims( dims_out, NULL);
	float Node_Generalized_Displacements[dims_out[0]][dims_out[1]]; Node_Generalized_Displacements_DataSet.read( Node_Generalized_Displacements, PredType::NATIVE_FLOAT);

	DataSet Node_Index_to_Generalized_Displacements_DataSet = file.openDataSet("Model/Nodes/Index_to_Generalized_Displacements");
	int Node_Index_to_Generalized_Displacements[Number_of_Nodes[0]+1]; Node_Index_to_Generalized_Displacements_DataSet.read( Node_Index_to_Generalized_Displacements, PredType::NATIVE_INT);

	DataSet Node_No_of_DOFs_DataSet = file.openDataSet("Model/Nodes/Number_of_DOFs");
	int Node_No_of_DOFs[Number_of_Nodes[0]+1]; Node_No_of_DOFs_DataSet.read( Node_No_of_DOFs, PredType::NATIVE_INT);


 	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// //////////////////////////////////////////////// Changed ////////////////////////////////////////
	// vtkDebugMacro( << "Send GMV data to Paraview");

	//  this->UpdateProgress(0.0);

	//    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
	//    {

	//    // Get the requested time step. We only support requests of a single time
	//    // step in this reader right now
	//    double requestedTimeValue =
	//      outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

	//    // Snapping and clamping of the requested time step to be in bounds is done by
	//    // FileSeriesReader class.
	//    vtkDebugMacro( << "RequestData: requested time value: " << requestedTimeValue);

	//    output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), requestedTimeValue);
	//    }

	// H5T_class_t type_class = dataset.getTypeClass();

	// if( type_class == H5T_INTEGER )
	// cout << "Data set has INTEGER type" << endl;
	// else 
	// cout << "Data set is different than Integers " << endl;

	// DataSpace dataspace = dataset.getSpace();
	// int rank = dataspace.getSimpleExtentNdims();
	// hsize_t dims_out[1];
	// int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
	// cout << "rank " << rank << ", dimensions " << (unsigned long)(dims_out[0])<< endl;

	// int nonode = dims_out[0];

	// /* memory space dimensions */
	// float  data_out[nonode]; hsize_t dimsm[1]; dimsm[0] = nonode; DataSpace memspace( 1, dimsm );

	// dataset.read( data_out, PredType::NATIVE_FLOAT, memspace, dataspace);

	vtkSmartPointer<vtkUnstructuredGrid> UGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkFloatArray> vtk_Generalized_Displacements = vtkSmartPointer<vtkFloatArray>::New();
	vtk_Generalized_Displacements->SetName("Generalized_Displacements");
	vtk_Generalized_Displacements->SetNumberOfComponents(3);
	vtk_Generalized_Displacements->SetComponentName(0,"X-axis");
	vtk_Generalized_Displacements->SetComponentName(1,"Y-axis");
	vtk_Generalized_Displacements->SetComponentName(2,"Z-axis");

	vtkSmartPointer<vtkFloatArray> vtk_Generalized_Velocity = vtkSmartPointer<vtkFloatArray>::New();
	vtk_Generalized_Velocity->SetName("Velocity");
	vtk_Generalized_Velocity->SetNumberOfComponents(3);
	vtk_Generalized_Velocity->SetComponentName(0,"X-axis");
	vtk_Generalized_Velocity->SetComponentName(1,"Y-axis");
	vtk_Generalized_Velocity->SetComponentName(2,"Z-axis");

	vtkSmartPointer<vtkFloatArray> vtk_Generalized_Acceleration = vtkSmartPointer<vtkFloatArray>::New();
	vtk_Generalized_Acceleration->SetName("Acceleration");
	vtk_Generalized_Acceleration->SetNumberOfComponents(3);
	vtk_Generalized_Acceleration->SetComponentName(0,"X-axis");
	vtk_Generalized_Acceleration->SetComponentName(1,"Y-axis");
	vtk_Generalized_Acceleration->SetComponentName(2,"Z-axis");

	vtkSmartPointer<vtkFloatArray> Elastic_Strain_Tensor = vtkSmartPointer<vtkFloatArray>::New();
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

	vtkSmartPointer<vtkFloatArray> Plastic_Strain_Tensor = vtkSmartPointer<vtkFloatArray>::New();
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

	vtkSmartPointer<vtkFloatArray> Stress_Tensor = vtkSmartPointer<vtkFloatArray>::New();
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

	vtkSmartPointer<vtkFloatArray> Material_Properties = vtkSmartPointer<vtkFloatArray>::New();
	Material_Properties->SetName("Material_Properties");
	Material_Properties->SetNumberOfComponents(1);

	vtkSmartPointer<vtkFloatArray> Total_Energy = vtkSmartPointer<vtkFloatArray>::New();
	Total_Energy->SetName("Total-Energy");
	Total_Energy->SetNumberOfComponents(1);

	vtkSmartPointer<vtkFloatArray> Incremental_Energy = vtkSmartPointer<vtkFloatArray>::New();
	Incremental_Energy->SetName("Incremental_Energy");
	Incremental_Energy->SetNumberOfComponents(1);


	points->SetNumberOfPoints(Number_of_Nodes[0]);
	UGrid->Allocate(Number_of_Elements[0]);
	float *pts  = (float *) points->GetVoidPointer(0);

	for (int i = 0; i <= Number_of_Nodes[0]; i++){
		if (Node_Index_to_Coordinates[i]!=-1){
			points->InsertPoint(i, Node_Coordinates[Node_Index_to_Coordinates[i]],Node_Coordinates[Node_Index_to_Coordinates[i]+1],Node_Coordinates[Node_Index_to_Coordinates[i]+2]);
			
			float tuple[3]={Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[i]][Ctime],Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[i]+1][Ctime],Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[i]+2][Ctime]};	
			vtk_Generalized_Displacements->InsertTupleValue(i,tuple);

			vtk_Generalized_Velocity->InsertTupleValue(i,tuple);
			vtk_Generalized_Acceleration->InsertTupleValue(i,tuple);
		}
	}

	UGrid->SetPoints(points);
	// cout << UGrid->GetNumberOfPoints() << endl;

	///////////////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	for (int i = 0; i <= Number_of_Elements[0]; i++){

		int connectivity_index = Element_Index_to_Connectivity[i];

		if(connectivity_index!=-1){
	
			int No_of_Element_Nodes = Element_Number_of_Nodes[i];
			vtkIdType Vertices[No_of_Element_Nodes];

			int Cell_Type = ESSI_to_VTK_Element.find(No_of_Element_Nodes)->second;
			std::vector<int> Nodes_Connectivity_Order = ESSI_to_VTK_Connectivity.find(No_of_Element_Nodes)->second;

			/************************************************** Gauss Calculations **********************************************/

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
				Elastic_Strain_Tensor->InsertNextTuple(El_Strain_Tuple);
				Plastic_Strain_Tensor->InsertNextTuple(Pl_Strain_Tuple);
				Stress_Tensor->InsertNextTuple(Stress_Tuple);
			}

			Incremental_Energy->InsertValue(i,this->Incremental_Energy_Database[Ctime*(Number_of_Elements[0]+1)+i]);
			Total_Energy->InsertValue(i,this->Total_Energy_Database[Ctime*(Number_of_Elements[0]+1)+i]);

			/********************************************************************************************************************/

			for(int j=0; j<No_of_Element_Nodes ; j++){
				Vertices[j] = Element_Connectivity[connectivity_index+Nodes_Connectivity_Order[j]];
			}

			UGrid->InsertNextCell(Cell_Type, No_of_Element_Nodes, Vertices);

			Material_Properties->InsertValue(i,i%10);
		}
	}
	
	/////////////////////////////////////////////////////////////////////// Building Up Data Array    ///////////////////////////////////////////////////////////////

	UGrid->GetPointData()->AddArray(vtk_Generalized_Displacements);
  	UGrid->GetPointData()->AddArray(vtk_Generalized_Velocity);
  	UGrid->GetPointData()->AddArray(vtk_Generalized_Acceleration);
 	UGrid->GetFieldData()->AddArray(Stress_Tensor);
  	UGrid->GetFieldData()->AddArray(Plastic_Strain_Tensor);
  	UGrid->GetFieldData()->AddArray(Elastic_Strain_Tensor);
  	UGrid->GetCellData()->AddArray(Material_Properties);
  	UGrid->GetCellData()->AddArray(Total_Energy);
  	UGrid->GetCellData()->AddArray(Incremental_Energy);

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

	output->ShallowCopy(UGrid);
	// UGrid->SetCells(NULL);
	// cout << UGrid->GetNumberOfCells() << endl;
	// cout << UGrid->GetNumberOfPoints() << endl;

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

