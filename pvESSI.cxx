#include "pvESSI.h"
 
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkCell.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkUnstructuredGrid.h"
#include "H5Cpp.h"
#include <string>
#include <iostream>
#include <map>
#include <vector>

// LINK_LIBRARIES(hdf5_cpp hdf5 )
// include_directories(/usr/include/hdf5/serial)

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
 
int pvESSI::RequestData(vtkInformation *vtkNotUsed(request),vtkInformationVector **vtkNotUsed(inputVector),	vtkInformationVector *outputVector){
 
	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

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

	// DataSet Element_Outputs_DataSet = file.openDataSet("Model/Elements/Outputs");
	// DataSpace Element_Outputs_DataSpace = Element_Outputs_DataSet.getSpace(); ndims = Element_Outputs_DataSpace.getSimpleExtentDims( dims_out, NULL);
	// float Element_Outputs[dims_out[0]][dims_out[1]]; Element_Outputs_DataSet.read( Element_Outputs, PredType::NATIVE_FLOAT);

	/////////////////////////////////////////////////// Reading Nodes datat //////////////////////////////////////////////////////////////////////////////////////////////

	DataSet Node_Coordinates_DataSet = file.openDataSet("Model/Nodes/Coordinates");
	float Node_Coordinates[Number_of_Nodes[0]*3]; Node_Coordinates_DataSet.read( Node_Coordinates, PredType::NATIVE_FLOAT);

	DataSet Node_Index_to_Coordinates_DataSet = file.openDataSet("Model/Nodes/Index_to_Coordinates");
	int Node_Index_to_Coordinates[Number_of_Nodes[0]+1]; Node_Index_to_Coordinates_DataSet.read( Node_Index_to_Coordinates, PredType::NATIVE_INT);

	// DataSet Node_Generalized_Displacements_DataSet = file.openDataSet("Model/Nodes/Generalized_Displacements");
	// DataSpace Node_Generalized_Displacements_DataSpace = Node_Generalized_Displacements_DataSet.getSpace(); ndims = Node_Generalized_Displacements_DataSpace.getSimpleExtentDims( dims_out, NULL);
	// float Node_Generalized_Displacements[dims_out[0]][dims_out[1]]; Node_Generalized_Displacements_DataSet.read( Node_Generalized_Displacements, PredType::NATIVE_FLOAT);

	DataSet Node_Index_to_Generalized_Displacements_DataSet = file.openDataSet("Model/Nodes/Index_to_Generalized_Displacements");
	int Node_Index_to_Generalized_Displacements[Number_of_Nodes[0]+1]; Node_Index_to_Generalized_Displacements_DataSet.read( Node_Index_to_Generalized_Displacements, PredType::NATIVE_INT);

	DataSet Node_No_of_DOFs_DataSet = file.openDataSet("Model/Nodes/Number_of_DOFs");
	int Node_No_of_DOFs[Number_of_Nodes[0]+1]; Node_No_of_DOFs_DataSet.read( Node_No_of_DOFs, PredType::NATIVE_INT);


 	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 	//  //  H5T_class_t type_class = dataset.getTypeClass();
	
	// // if( type_class == H5T_INTEGER )
	// // 	cout << "Data set has INTEGER type" << endl;
	// // else 
	// // 	cout << "Data set is different than Integers " << endl;

	// // DataSpace dataspace = dataset.getSpace();
	// // int rank = dataspace.getSimpleExtentNdims();
	// // hsize_t dims_out[1];
	// // int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
	// // cout << "rank " << rank << ", dimensions " << (unsigned long)(dims_out[0])<< endl;

	// // int nonode = dims_out[0];

 	//  //  /* memory space dimensions */
	// // float  data_out[nonode]; hsize_t dimsm[1]; dimsm[0] = nonode; DataSpace memspace( 1, dimsm );

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

	///////////////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	for (int i = 0; i <= Number_of_Elements[0]; i++){

		int connectivity_index = Element_Index_to_Connectivity[i];

		if(connectivity_index==-1)
			continue;

		int No_of_Element_Nodes = Element_Number_of_Nodes[i];

		vtkIdType Vertices[No_of_Element_Nodes];

		int Cell_Type = ESSI_to_VTK_Element.find(No_of_Element_Nodes)->second;
		std::vector<int> Nodes_Connectivity_Order = ESSI_to_VTK_Connectivity.find(No_of_Element_Nodes)->second;

		for(int j=0; j<No_of_Element_Nodes ; j++)
			Vertices[j] = Element_Connectivity[connectivity_index+Nodes_Connectivity_Order[j]];

		UGrid->InsertNextCell(Cell_Type, No_of_Element_Nodes, Vertices);
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





	/****************************************************************************************************/

	// Here is where you would read the data from the file. In this example,
	// we simply create a point.

	
	// polydata->SetPoints(points);

	// vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
	// glyphFilter->SetInputData(UGrid);
	// glyphFilter->Update();

	output->ShallowCopy(UGrid);
	// cout << UGrid->GetNumberOfCells() << endl;
	// cout << UGrid->GetNumberOfPoints() << endl;

	return 1;
}
 
void pvESSI::PrintSelf(ostream& os, vtkIndent indent){

	this->Superclass::PrintSelf(os,indent);
	os << indent << "File Name: " << (this->FileName ? this->FileName : "(none)") << "\n";
}

void pvESSI::set_VTK_To_ESSI_Elements_Connectivity(){

	std::vector<int> connectivity_vector;

	ESSI_to_VTK_Element[1] = VTK_VERTEX;		ESSI_to_VTK_Element[2]  = VTK_LINE;		ESSI_to_VTK_Element[4] = VTK_QUAD;
	ESSI_to_VTK_Element[8] = VTK_HEXAHEDRON;	/*ESSI_to_VTK_Element[20] = VTK_VERTEX;*/	ESSI_to_VTK_Element[27] = VTK_TRIQUADRATIC_HEXAHEDRON;

	int Node[1] = {0};						connectivity_vector.assign(Node,Node+1);				ESSI_to_VTK_Connectivity[1] = connectivity_vector;	connectivity_vector.clear();
	int Line[2] = {0,1}; 					connectivity_vector.assign(Line,Line+2); 				ESSI_to_VTK_Connectivity[2] = connectivity_vector;	connectivity_vector.clear();
	int Quadrangle[4] = {0,1,2,3}; 			connectivity_vector.assign(Quadrangle,Quadrangle+4); 	ESSI_to_VTK_Connectivity[4] = connectivity_vector;	connectivity_vector.clear();
	int Hexahedron[8] = {4,5,6,7,0,1,2,3}; 	connectivity_vector.assign(Hexahedron,Hexahedron+8);	ESSI_to_VTK_Connectivity[8] = connectivity_vector;	connectivity_vector.clear();
	int TriQudratic_hexahedron[27] = {6,5,4,7,2,1,0,3,13,12,15,14,9,8,11,10,18,17,16,19,23,21,22,24,26,25,20}; 											connectivity_vector.assign(TriQudratic_hexahedron,TriQudratic_hexahedron+27);					ESSI_to_VTK_Connectivity[27]= connectivity_vector;	connectivity_vector.clear();
}