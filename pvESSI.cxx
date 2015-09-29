#include "pvESSI.h"
 
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
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
}
 


// int pvESSI::RequestData(vtkInformation *request,vtkInformationVector **vtkNotUsed(inputVector),	vtkInformationVector *outputVector){
 

//   if(request->Get( vtkDemandDrivenPipeline::FROM_OUTPUT_PORT() ) < 0)
//     {
//     this->GetExecutive()->GetOutputData(0)->Initialize();
//     return 0;
//     }

//   this->CurrentTimeStep = this->TimeStep;

//   // Get the output pipeline information and data object.
//   vtkInformation* outInfo = outputVector->GetInformationObject(0);
//   vtkDataObject* output = outInfo->Get(vtkDataObject::DATA_OBJECT());

//   // Check if a particular time was requested.
//   if(outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
//     {
//     // Get the requested time step.
//     this->CurrentTimeStep =
//       outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

//     // Save the time value in the output data information.
//     int length =
//       outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
//     if(this->CurrentTimeStep >= 0 && this->CurrentTimeStep < length)
//       {
//       double* steps =
//         outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
//       // output->GetInformation()->Set(vtkDataObject::DATA_TIME(),
//       //                               steps[this->CurrentTimeStep]);
//       }
//     else
//       {
//       vtkErrorMacro("Time index " << this->CurrentTimeStep
//                     << " requested but there are "
//                     << length << " time steps.");
//       }

//     // Clamp the requested time step to be in bounds.
//     if ( this->CurrentTimeStep < this->TimeStepRange[0] )
//       {
//       this->CurrentTimeStep = this->TimeStepRange[0];
//       }
//     else if ( this->CurrentTimeStep > this->TimeStepRange[1] )
//       {
//       this->CurrentTimeStep = this->TimeStepRange[1];
//       }
//     }

//   // Set the time we will store in the output.
//   // output->GetInformation()->Set(vtkDataObject::DATA_TIME_INDEX(),
//   //                               this->CurrentTimeStep);

//   // Re-open the input file.  If it fails, the error was already
//   // reported by OpenVTKFile.

// return 1;

// }
int pvESSI::RequestData(vtkInformation *vtkNotUsed(request),vtkInformationVector **vtkNotUsed(inputVector),	vtkInformationVector *outputVector){
 
	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	outInfo->Print(std::cout);

	int extent[6] = {0,-1,0,-1,0,-1};
  	outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);

  	double Ctime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
  	int Clength = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  	double* Csteps = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

  	cout<< "Ctime " << "  " << Ctime << endl;
  	cout<< "Clength " << "  " << Clength << endl;
  	for (int i =0 ; i< Clength ; i++){
  		cout << Csteps[i] << "  " << endl;
  	}

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




 		////////////////////////////////////////////////// Changed ////////////////////////////////////////
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

	points->SetNumberOfPoints(Number_of_Nodes[0]);
	UGrid->Allocate(Number_of_Elements[0]);
	float *pts  = (float *) points->GetVoidPointer(0);

	for (int i = 0; i <= Number_of_Nodes[0]; i++){
		if (Node_Index_to_Coordinates[i]>-1){
			points->InsertPoint(i, Node_Coordinates[Node_Index_to_Coordinates[i]],Node_Coordinates[Node_Index_to_Coordinates[i]+1],Node_Coordinates[Node_Index_to_Coordinates[i]+2]);
		}
	}

	UGrid->SetPoints(points);
	// cout << UGrid->GetNumberOfPoints() << endl;

	this->set_VTK_To_ESSI_Elements_Connectivity();

		vtkSmartPointer<vtkFloatArray> vtk_Generalized_Displacements = vtkSmartPointer<vtkFloatArray>::New();
	vtk_Generalized_Displacements->SetName("Generalized_Displacements");
	vtk_Generalized_Displacements->SetNumberOfComponents(3);

	///////////////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	for (int i = 0; i <= Number_of_Elements[0]; i++){

		int connectivity_index = Element_Index_to_Connectivity[i];

		if(connectivity_index==-1)
			continue;

		float tuple[3];
		tuple[0]=i;
		tuple[1]=i+1;
		tuple[2]=i+2;

		int No_of_Element_Nodes = Element_Number_of_Nodes[i];
		vtk_Generalized_Displacements->InsertNextTuple(tuple);
		vtk_Generalized_Displacements->SetComponentName(0,"X-axis");
		vtk_Generalized_Displacements->SetComponentName(1,"Y-axis");
		vtk_Generalized_Displacements->SetComponentName(2,"Z-axis");

		vtkIdType Vertices[No_of_Element_Nodes];

		int Cell_Type = ESSI_to_VTK_Element.find(No_of_Element_Nodes)->second;
		std::vector<int> Nodes_Connectivity_Order = ESSI_to_VTK_Connectivity.find(No_of_Element_Nodes)->second;

		for(int j=0; j<No_of_Element_Nodes ; j++)
			Vertices[j] = Element_Connectivity[connectivity_index+Nodes_Connectivity_Order[j]];

		UGrid->InsertNextCell(Cell_Type, No_of_Element_Nodes, Vertices);
	}
	
	/////////////////////////////////////////////////////////////////////// Building Up Data Array    ///////////////////////////////////////////////////////////////

	UGrid->GetCellData()->AddArray(vtk_Generalized_Displacements);


	/****************************************************************************************************/

	// Here is where you would read the data from the file. In this example,
	// we simply create a point.

	
	// polydata->SetPoints(points);

	// vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
	// glyphFilter->SetInputData(UGrid);
	// glyphFilter->Update();

	output->ShallowCopy(UGrid);
	// UGrid->SetCells(NULL);
	// cout << UGrid->GetNumberOfCells() << endl;
	// cout << UGrid->GetNumberOfPoints() << endl;

	return 1;
}

int pvESSI::RequestInformation( vtkInformation *request, vtkInformationVector **vtkNotUsed(inVec), vtkInformationVector* outVec){

	vtkInformation* outInfo = outVec->GetInformationObject(0);
	int extent[6]={0,100,0,100,0,100}; 
	int  No_of_TimeSteps = 20;
	this->NumberOfTimeSteps = No_of_TimeSteps;
	double Time_Steps[No_of_TimeSteps];

	for (int i =0;i<No_of_TimeSteps ;i++)
		Time_Steps[i] = i;
	double Time_range[2]={0,No_of_TimeSteps-1};
	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),extent,6);
	outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),Time_Steps,No_of_TimeSteps);
	outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),Time_range,2);
	// outInfo->Set(vtkAlgorithm::CAN_PRODUCE_SUB_EXTENT(),1);


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