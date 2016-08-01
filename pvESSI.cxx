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

// LINK_LIBRARIES(hdf5_cpp hdf5 )

/************************************************************************************************************************************************/
// cmake .. -DParaView_DIR=~/Softwares/Paraview/Paraview-Build/ -DGHOST_BUILD_CDAWEB=OFF
// cmake .. -DParaView_DIR="/home/sumeet/Softwares/ParaView-v4.4.0" -DGHOST_BUILD_CDAWEB=OFF
/************************************************************************************************************************************************/


vtkStandardNewMacro(pvESSI);

pvESSI::pvESSI(){ 


	this->FileName = NULL;
	this->SetNumberOfInputPorts(0);
	this->SetNumberOfOutputPorts(2);
	UGrid_Gauss_Mesh          = vtkSmartPointer<vtkUnstructuredGrid>::New();
	UGrid_Node_Mesh           = vtkSmartPointer<vtkUnstructuredGrid>::New();
	UGrid_Current_Node_Mesh   = vtkSmartPointer<vtkUnstructuredGrid>::New();
	UGrid_Current_Gauss_Mesh  = vtkSmartPointer<vtkUnstructuredGrid>::New();
	this->set_VTK_To_ESSI_Elements_Connectivity();
	this->Build_Meta_Array_Map();
	this->initialize();
}

/*****************************************************************************
* This method responds to the request made by vtk Pipeleine. 
* This method is invoked when the time stamp is changed from paraview VCR.
*****************************************************************************/

int pvESSI::RequestData(vtkInformation *vtkNotUsed(request),vtkInformationVector **vtkNotUsed(inputVector),	vtkInformationVector *outputVector){
 
	// get the info object
	vtkInformation *Node_Mesh = outputVector->GetInformationObject(0);
	vtkInformation *Gauss_Mesh = outputVector->GetInformationObject(1);
	// outInfo->Print(std::cout);

	if (!Whether_Node_Mesh_Build){
		this->Get_Node_Mesh(UGrid_Node_Mesh);
		UGrid_Current_Node_Mesh->ShallowCopy(UGrid_Node_Mesh);
	} 

	// this->Get_Node_Mesh();
	// this->Get_Gauss_Mesh();
	// 
	if (!Whether_Gauss_Mesh_Build ){ 
		this->Get_Gauss_Mesh(UGrid_Gauss_Mesh);
		UGrid_Current_Gauss_Mesh->ShallowCopy(UGrid_Gauss_Mesh);
	} 

	// int extent[6] = {0,-1,0,-1,0,-1};
	// outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);

  	this->Node_Mesh_Current_Time = Time_Map.find( Node_Mesh->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))->second;
  	// Gauss_Mesh_Current_Time = Node_Mesh_Current_Time;
  	this->Gauss_Mesh_Current_Time = Time_Map.find( Gauss_Mesh->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))->second;

  	// int Clength = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  	// double* Csteps = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

  	if(Gauss_Mesh_Current_Time > Number_Of_Time_Steps){
  		Gauss_Mesh_Current_Time =0;
  	}


 	std::cout<< "Node_Mesh_Current_Time " << Node_Mesh_Current_Time<< std::endl;
 	std::cout<< "Gauss_Mesh_Current_Time " << Gauss_Mesh_Current_Time<< std::endl;

 	// vtkIndent indent;
 	// outputVector->PrintSelf(std::cout, indent);
 	// delaunay3D->PrintSelf(std::cout, indent);

 // 	cout << "Number of Nodes "   << " " << Number_of_Nodes << endl;
	// cout << "Pseudo_Number_of_Nodes "  << " " << Pseudo_Number_of_Nodes << endl;
	// cout << " Number of Gauss Points " << Number_of_Gauss_Nodes << endl;

  // Generate a tetrahedral mesh from the input points. By
  // default, the generated volume is the convex hull of the points.

  //Read the file



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
		Build_Node_Attributes(UGrid_Current_Node_Mesh, this->Node_Mesh_Current_Time );
	// 	Node_Mesh_Previous_Time = Node_Mesh_Current_Time;
	// }
	// else if(Gauss_Mesh_Current_Time!=Gauss_Mesh_Previous_Time){
		Build_Gauss_Attributes(UGrid_Current_Gauss_Mesh, this->Gauss_Mesh_Current_Time );
	// 	Gauss_Mesh_Previous_Time = Gauss_Mesh_Current_Time;
	// }

	vtkUnstructuredGrid *Output_Node_Mesh = vtkUnstructuredGrid::SafeDownCast(Node_Mesh->Get(vtkDataObject::DATA_OBJECT()));
	vtkUnstructuredGrid *Output_Gauss_Mesh = vtkUnstructuredGrid::SafeDownCast(Gauss_Mesh->Get(vtkDataObject::DATA_OBJECT()));
	Output_Node_Mesh->ShallowCopy(UGrid_Current_Node_Mesh);
	Output_Gauss_Mesh->ShallowCopy(UGrid_Current_Gauss_Mesh);

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
	// vtkPolyData *Coutput=ildvtkPolyData::SafeDownCast(pointGen->GetOutput());
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

/*****************************************************************************
* This method is called only once for information about time stamp and extent. 
* This is the method which is called firt after the default constructor is 
* initialized
*****************************************************************************/

int pvESSI::RequestInformation( vtkInformation *request, vtkInformationVector **vtkNotUsed(inVec), vtkInformationVector* outVec){

	vtkInformation* Node_Mesh = outVec->GetInformationObject(0);
	vtkInformation* Gauss_Mesh = outVec->GetInformationObject(1);

	herr_t status;
	hid_t File, DataSet, DataSpace; 
	hsize_t  dims_out[1];

	//////////////////////////////////////////////////////// Reading Hdf5 File ///////////////////////////////////////////////////////////////////////////////////////

	File = H5Fopen(this->FileName, H5F_ACC_RDWR, H5P_DEFAULT);
	int Build_Map_Status[1]={1};

	////////////////////////////////////////////////////// Reading General Outside Data /////////////////////////////////////////////////////////////////////////////

	DataSet = H5Dopen(File, "Number_of_Time_Steps", H5P_DEFAULT);
	int No_of_TimeSteps[1]; status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,No_of_TimeSteps);
	this->Number_Of_Time_Steps = No_of_TimeSteps[0]; status = H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Number_of_Elements", H5P_DEFAULT);
	int Number_of_Elements[1]; status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Number_of_Elements);
	this->Number_of_Elements = Number_of_Elements[0]; status = H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Number_of_Nodes", H5P_DEFAULT);
	int Number_of_Nodes[1]; status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Number_of_Nodes);
	this->Number_of_Nodes = Number_of_Nodes[0]; status = H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Whether_Maps_Build", H5P_DEFAULT);

	if(DataSet < 0){
		/*** Dataset not present so lets create it **/
		dims_out[0]=1; DataSpace = H5Screate_simple(1, dims_out, NULL);
		DataSet = H5Dcreate(File,"Whether_Maps_Build",H5T_STD_I32BE,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
		Build_Map_Status[0] = -1; status = H5Sclose(DataSpace); status = H5Dclose(DataSet); 
	}
	else{
		status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Build_Map_Status);
	    status = H5Dclose(DataSet);
	}

   	if(Build_Map_Status[0] <= 0){
   		
   		Build_Maps();	/*** Maps Not Built So Lets Go and Create the Map Dataset **/

   	}
   	else{

   		DataSet = H5Dopen(File, "Model/Nodes/Index_to_Coordinates", H5P_DEFAULT); DataSpace = H5Dget_space(DataSet);
		status  = H5Sget_simple_extent_dims(DataSpace, dims_out, NULL);this->Pseudo_Number_of_Nodes = dims_out[0];
		status = H5Sclose(DataSpace); status = H5Dclose(DataSet);

   		DataSet = H5Dopen(File, "Model/Elements/Index_to_Connectivity", H5P_DEFAULT); DataSpace = H5Dget_space(DataSet);
		status  = H5Sget_simple_extent_dims(DataSpace, dims_out, NULL);this->Pseudo_Number_of_Elements = dims_out[0];
		status = H5Sclose(DataSpace); status = H5Dclose(DataSet);
   }

	///////////////////////////////////////////////////////////////// Evaluating Number of Gauss Points ////////////////////////////////////////////////////////////////////////////////////////////

   	DataSet = H5Dopen(File, "Model/Elements/Number_of_Gauss_Points", H5P_DEFAULT); int Element_Number_of_Gauss_Points[Pseudo_Number_of_Elements]; 
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Number_of_Gauss_Points); status = H5Dclose(DataSet);
	this->Number_of_Gauss_Nodes = 0;

    for (int i=0; i < this->Pseudo_Number_of_Elements; i++){
        if (Element_Number_of_Gauss_Points[i] > 0){
            this->Number_of_Gauss_Nodes += Element_Number_of_Gauss_Points[i];
        }
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	DataSet = H5Dopen(File, "time", H5P_DEFAULT); this->Time = new double[Number_Of_Time_Steps]; 
	status = H5Dread(DataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,Time);  status = H5Dclose(DataSet);

	Build_Time_Map(); 

	double Time_range[2]={Time[0],Time[Number_Of_Time_Steps-1]};

	
	Node_Mesh->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),Time, No_of_TimeSteps[0]);
	Node_Mesh->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),Time_range,2);

	Gauss_Mesh->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),Time, No_of_TimeSteps[0]);
	Gauss_Mesh->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),Time_range,2);

	// outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),extent,6);
	// outInfo->Set(vtkAlgorithm::CAN_PRODUCE_SUB_EXTENT(),1);


	/////////////////////////////////////////////// Assign Key /////////////////////////////////////////////////////
	// vtkSmartPointer<vtkInformationQuadratureSchemeDefinitionVectorKey> Quadrature_Definition_Vector = vtkSmartPointer<vtkInformationQuadratureSchemeDefinitionVectorKey> ::New();

	// vtkSmartPointer<vtkQuadratureSchemeDefinition> Quadrature_Definition = vtkSmartPointer<vtkQuadratureSchemeDefinition>::New();
	// Quadrature_Definition->Initialize(VTK_HEXAHEDRON,8,8,W_QQ_64_A);

	status = H5Fclose(File);
	return 1;
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

	herr_t status;
	hid_t File, DataSet, DataSpace, MemSpace; 
	hsize_t  dims_out[2], offset[2], count[2], col_dims[1];

	//////////////////////////////////////////////////////// Reading Hdf5 File ///////////////////////////////////////////////////////////////////////////////////////

	File = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);

	//////////////////////////////////////////////////// Reading Map Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet = H5Dopen(File, "Maps/NodeMap", H5P_DEFAULT); int *Node_Map; Node_Map = (int*) malloc((Number_of_Nodes) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Node_Map); status=H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Maps/ElementMap", H5P_DEFAULT); int *Element_Map; Element_Map = (int*) malloc((Number_of_Elements) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Map); status=H5Dclose(DataSet);

	////////////////////////////////////////////////////// Reading Node Attributes /////////////////////////////////////////////////////////////////////////////

	DataSet = H5Dopen(File, "Model/Nodes/Index_to_Generalized_Displacements", H5P_DEFAULT); int *Node_Index_to_Generalized_Displacements; Node_Index_to_Generalized_Displacements = (int*) malloc((Pseudo_Number_of_Nodes) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Node_Index_to_Generalized_Displacements); status=H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Model/Elements/Material_tags", H5P_DEFAULT); int *Element_Material_Tags; Element_Material_Tags = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Material_Tags); status=H5Dclose(DataSet);

	///////////////////////////////////////////  Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	

	DataSet = H5Dopen(File, "Model/Nodes/Generalized_Displacements", H5P_DEFAULT); DataSpace = H5Dget_space(DataSet);
	status  = H5Sget_simple_extent_dims(DataSpace, dims_out, NULL);	float *Node_Generalized_Displacements; Node_Generalized_Displacements = (float*) malloc((dims_out[0] ) * sizeof(float));
	offset[0]=0; 					  count[0] = dims_out[0];		col_dims[0]=dims_out[0];
	offset[1]=Current_Time; 		  count[1] = 1;					MemSpace = H5Screate_simple(1,col_dims,NULL);
	status = H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset,NULL,count,NULL);
	status = H5Dread(DataSet, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, Node_Generalized_Displacements); 
	status = H5Sclose(MemSpace); status=H5Sclose(DataSpace); status=H5Dclose(DataSet);

 	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Generalized_Displacements = vtkSmartPointer<vtkFloatArray>::New(); 
	this->Set_Meta_Array (Meta_Array_Map["Generalized_Displacements"]);

	Material_Tag = vtkSmartPointer<vtkIntArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Material_Tag"]);

	Node_Tag = vtkSmartPointer<vtkIntArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Node_Tag"]);

	Element_Tag = vtkSmartPointer<vtkIntArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Element_Tag"]);	

 	/////////////////////////////////////////////////////////////// DataSets Visulization at Nodes //////////////////////////////////////////////////////////////////////////////////////
 	int Node_No; 
	for (int i = 0; i < Number_of_Nodes; i++){

			Node_No = Node_Map[i]; Node_Tag -> InsertValue(i,Node_No);

			float tuple[3]={
				Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[Node_No]  ],
				Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[Node_No]+1],
				Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[Node_No]+2]
			};	

			Generalized_Displacements->InsertTypedTuple(i,tuple);


			// ************************************************************** Displacement, Acceleration and Velocity -> Calculation Formulae *********************************************************

			// Displacement(t) = D_t;
			// Acceleration(t) = (1/del.t)^2*{ D_(t+del.t) - 2*D_t + D_(t-del.t)};
			// Acceleration(t) = 1/2/del.t*{ D_(t+del.t) - D_(t-del.t)};

			/***************************************************************** Add Acceleration and Velocity Arrays Here ****************************************/
				
			/****************************************************************************************************************************************************/

	}

	Node_Mesh->GetPointData()->AddArray(Generalized_Displacements);
	Node_Mesh->GetPointData()->AddArray(Node_Tag);

	// /////////////////////////////////////////////////////////////////////// DataSet Visualization in Cell   ///////////////////////////////////////////////////////////////
 	int Element_No; 
	for (int i = 0; i < Number_of_Elements; i++){

		Element_No = Element_Map[i]; 
		Element_Tag -> InsertValue(i,Element_No);
		Material_Tag -> InsertValue(i,Element_Material_Tags[Element_No]);

	}

	Node_Mesh->GetCellData()->AddArray(Element_Tag);
	Node_Mesh->GetCellData()->AddArray(Material_Tag);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	H5Fclose(File);

	// freeing heap allocation for variables 
	free(Node_Map);
  	free(Element_Map);
  	free(Node_Index_to_Generalized_Displacements);
  	free(Element_Material_Tags);
	free(Node_Generalized_Displacements);
	
  	return;
}

/*****************************************************************************
* This Method builds gauss attributes for the current time and pushes it to 
* the <vtkUnstructuredGrid> input vtkobject.
*****************************************************************************/

void pvESSI::Build_Gauss_Attributes(vtkSmartPointer<vtkUnstructuredGrid> Gauss_Mesh, int Current_Time){

	this->Build_ProbeFilter_Gauss_Mesh(Gauss_Mesh, 0);

	herr_t status;
	hid_t File, DataSet, DataSpace, MemSpace; 
	hsize_t  dims_out[2], offset[2], count[2], col_dims[1];

	//////////////////////////////////////////////////////// Reading Hdf5 File ///////////////////////////////////////////////////////////////////////////////////////

	File = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);

	//////////////////////////////////////////////////// Reading Map Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet = H5Dopen(File, "Maps/ElementMap", H5P_DEFAULT); int *Element_Map; Element_Map = (int*) malloc((Number_of_Elements) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Map); status=H5Dclose(DataSet);

	////////////////////////////////////////////////////// Reading Element  Attributes /////////////////////////////////////////////////////////////////////////////

	DataSet = H5Dopen(File, "Model/Elements/Index_to_Outputs", H5P_DEFAULT); int *Element_Index_to_Outputs; Element_Index_to_Outputs = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Index_to_Outputs); status=H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Model/Elements/Number_of_Gauss_Points", H5P_DEFAULT); int *Element_Number_of_Gauss_Points; Element_Number_of_Gauss_Points = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Number_of_Gauss_Points); status=H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Model/Elements/Number_of_Output_Fields", H5P_DEFAULT); int *Element_Number_of_Output_Fields; Element_Number_of_Output_Fields = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Number_of_Output_Fields); status=H5Dclose(DataSet);


	///////////////////////////////////////////  Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	

	DataSet = H5Dopen(File, "Model/Elements/Outputs", H5P_DEFAULT); DataSpace = H5Dget_space(DataSet);
	status  = H5Sget_simple_extent_dims(DataSpace, dims_out, NULL);	float *Element_Outputs; Element_Outputs = (float*) malloc((dims_out[0] ) * sizeof(float));
	offset[0]=0; 					   count[0] = dims_out[0];		col_dims[0]=dims_out[0];
	offset[1]=Current_Time; 		   count[1] = 1;					MemSpace = H5Screate_simple(1,col_dims,NULL);
	status = H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset,NULL,count,NULL);
	status = H5Dread(DataSet, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, Element_Outputs); 
	status = H5Sclose(MemSpace); status=H5Sclose(DataSpace); status=H5Dclose(DataSet);

	///////////////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	Elastic_Strain = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Elastic_Strain"]);

	Plastic_Strain = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Plastic_Strain"]);

	Stress = vtkSmartPointer<vtkFloatArray> ::New();
	this->Set_Meta_Array (Meta_Array_Map["Stress"]);


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

 	int Index_of_Gauss_Nodes =0, element_no, No_of_Element_Gauss_Nodes,Output_Index, gauss_number=0;

	for (int i = 0; i < Number_of_Elements; i++){

		element_no = Element_Map[i];

		No_of_Element_Gauss_Nodes = Element_Number_of_Gauss_Points[element_no];
		Output_Index = Element_Index_to_Outputs[element_no];

		for(int j=0; j< No_of_Element_Gauss_Nodes ; j++){

			float El_Strain_Tuple[6] ={
				Element_Outputs[Output_Index  ],
				Element_Outputs[Output_Index+3],
				Element_Outputs[Output_Index+4],
				Element_Outputs[Output_Index+1],
				Element_Outputs[Output_Index+5],
				Element_Outputs[Output_Index+2]
			};

			Output_Index = Output_Index+6;

			float Pl_Strain_Tuple[6] ={
				Element_Outputs[Output_Index  ],
				Element_Outputs[Output_Index+3],
				Element_Outputs[Output_Index+4],
				Element_Outputs[Output_Index+1],
				Element_Outputs[Output_Index+5],
				Element_Outputs[Output_Index+2]
			};

			Output_Index = Output_Index+6;

			float Stress_Tuple[6] ={
				Element_Outputs[Output_Index  ],
				Element_Outputs[Output_Index+3],
				Element_Outputs[Output_Index+4],
				Element_Outputs[Output_Index+1],
				Element_Outputs[Output_Index+5],
				Element_Outputs[Output_Index+2]
			};

			Output_Index = Output_Index+6;

			Elastic_Strain->InsertTypedTuple (Index_of_Gauss_Nodes,El_Strain_Tuple);
			Plastic_Strain->InsertTypedTuple (Index_of_Gauss_Nodes,Pl_Strain_Tuple);
			Stress->InsertTypedTuple (Index_of_Gauss_Nodes,Stress_Tuple);
			
			Index_of_Gauss_Nodes+=1;

		}
	}

	Gauss_Mesh->GetPointData()->AddArray(Stress);
  	Gauss_Mesh->GetPointData()->AddArray(Plastic_Strain);
  	Gauss_Mesh->GetPointData()->AddArray(Elastic_Strain);

	H5Fclose(File);

	// freeing heap allocation for variables 
  	free(Element_Map);
  	free(Element_Index_to_Outputs);
  	free(Element_Number_of_Gauss_Points);
  	free(Element_Number_of_Output_Fields);
	free(Element_Outputs);
	

	return;
}

/*****************************************************************************
* This Method probes node mesh variables at gauss mesh 
* By Default [prob_type =0], probes only generalized displacement
* But the user can choose an option [prob_type = 1] to probe all variables
*****************************************************************************/

void pvESSI::Build_ProbeFilter_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Probe_Input, int probe_type){

	vtkSmartPointer<vtkUnstructuredGrid> Probe_Source = vtkSmartPointer<vtkUnstructuredGrid>::New();
	Probe_Source->ShallowCopy(this->UGrid_Node_Mesh);

	if (probe_type){
		Build_Node_Attributes(Probe_Source, this->Gauss_Mesh_Current_Time);
	}
	else{

		herr_t status;
		hid_t File, DataSet, DataSpace, MemSpace; 
		hsize_t  dims_out[2], offset[2], count[2], col_dims[1];

		//////////////////////////////////////////////////////// Reading Hdf5 File ///////////////////////////////////////////////////////////////////////////////////////

		File = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);

		//////////////////////////////////////////////////// Reading Map Data ///////////////////////////////////////////////////////////////////////////////////////////

		DataSet = H5Dopen(File, "Maps/NodeMap", H5P_DEFAULT); int *Node_Map; Node_Map = (int*) malloc((Number_of_Nodes) * sizeof(int));
		status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Node_Map); status=H5Dclose(DataSet);

		////////////////////////////////////////////////////// Reading Node Attributes /////////////////////////////////////////////////////////////////////////////

		DataSet = H5Dopen(File, "Model/Nodes/Index_to_Generalized_Displacements", H5P_DEFAULT); int *Node_Index_to_Generalized_Displacements; Node_Index_to_Generalized_Displacements = (int*) malloc((Pseudo_Number_of_Nodes) * sizeof(int));
		status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Node_Index_to_Generalized_Displacements); status=H5Dclose(DataSet);

		///////////////////////////////////////////  Output Dataset for a particular time /////////////////////////////////////////////////////////////////////////////	

		DataSet = H5Dopen(File, "Model/Nodes/Generalized_Displacements", H5P_DEFAULT); DataSpace = H5Dget_space(DataSet);
		status  = H5Sget_simple_extent_dims(DataSpace, dims_out, NULL);	float *Node_Generalized_Displacements; Node_Generalized_Displacements = (float*) malloc((dims_out[0] ) * sizeof(float));
		offset[0]=0; 					  		 count[0] = dims_out[0];		col_dims[0]=dims_out[0];
		offset[1]=this->Gauss_Mesh_Current_Time; count[1] = 1;					MemSpace = H5Screate_simple(1,col_dims,NULL);
		status = H5Sselect_hyperslab(DataSpace,H5S_SELECT_SET,offset,NULL,count,NULL);
		status = H5Dread(DataSet, H5T_NATIVE_FLOAT, MemSpace, DataSpace, H5P_DEFAULT, Node_Generalized_Displacements); 
		status = H5Sclose(MemSpace); status=H5Sclose(DataSpace); status=H5Dclose(DataSet);

	 	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		Generalized_Displacements = vtkSmartPointer<vtkFloatArray>::New(); 
		this->Set_Meta_Array (Meta_Array_Map["Generalized_Displacements"]);
 	
	 	/////////////////////////////////////////////////////////////// DataSets Visulization at Nodes //////////////////////////////////////////////////////////////////////////////////////
	 	int Node_No; 
		for (int i = 0; i < Number_of_Nodes; i++){

				Node_No = Node_Map[i]; 
				float tuple[3]={
					Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[Node_No]],
					Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[Node_No]+1],
					Node_Generalized_Displacements[Node_Index_to_Generalized_Displacements[Node_No]+2]
				};	

				Generalized_Displacements->InsertTypedTuple(i,tuple);


		}

		Probe_Source->GetPointData()->AddArray(Generalized_Displacements);

		H5Fclose(File);

		// freeing heap allocation for variables  
		free(Node_Map);
	  	free(Node_Index_to_Generalized_Displacements);
		free(Node_Generalized_Displacements);

	}

	/************* Initializing Probe filter ******************************************/

	vtkSmartPointer<vtkProbeFilter> ProbeFilter = vtkSmartPointer<vtkProbeFilter>::New();
	ProbeFilter->SetSourceData(Probe_Source);
	ProbeFilter->SetInputData(Probe_Input);
	ProbeFilter->Update();

	Probe_Input->ShallowCopy( ProbeFilter->GetOutput());

	////////////////////////////////////// For Debugging ////////////////////////////////
	// UGrid_Gauss_Mesh->GetPointData()->RemoveArray("Material_Tag");
	// UGrid_Gauss_Mesh->GetPointData()->RemoveArray("Node_No");
	// UGrid_Gauss_Mesh->GetPointData()->RemoveArray("Element_No");
	/////////////////////////////////////////////////////////////////////////////////////

  	return;

}


/*****************************************************************************
* This function creates a Node Mesh i.e a mesh from the given node data 
* and element data from hdf5 file. The function stores the mesh in the 
* given input <vtkUnstructuredGrid> input object.
*****************************************************************************/

void pvESSI::Get_Node_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Node_Mesh){

	herr_t status;
	hid_t File, DataSet, DataSpace; 
	hsize_t  dims_out[1];

	//////////////////////////////////////////////////////// Reading Hdf5 File ///////////////////////////////////////////////////////////////////////////////////////

	File = H5Fopen(this->FileName, H5F_ACC_RDWR, H5P_DEFAULT);

	//////////////////////////////////////////////////// Reading Map Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet = H5Dopen(File, "Maps/NodeMap", H5P_DEFAULT); int *Node_Map; Node_Map = (int*) malloc((Number_of_Nodes) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Node_Map); status=H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Maps/InverseNodeMap", H5P_DEFAULT); int *Inverse_Node_Map; Inverse_Node_Map = (int*) malloc((Pseudo_Number_of_Nodes) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Inverse_Node_Map); status=H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Maps/ElementMap", H5P_DEFAULT); int *Element_Map; Element_Map = (int*) malloc((Number_of_Elements) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Map); status=H5Dclose(DataSet);

	//////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet = H5Dopen(File, "Model/Elements/Index_to_Connectivity", H5P_DEFAULT); int *Element_Index_to_Connectivity; Element_Index_to_Connectivity = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Index_to_Connectivity); status=H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Model/Elements/Class_Tags", H5P_DEFAULT); int *Element_Class_Tags; Element_Class_Tags = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Class_Tags); status=H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Model/Elements/Number_of_Nodes", H5P_DEFAULT); int *Element_Number_of_Nodes; Element_Number_of_Nodes = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Number_of_Nodes); status=H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Model/Elements/Connectivity", H5P_DEFAULT); DataSpace = H5Dget_space(DataSet);
	status  = H5Sget_simple_extent_dims(DataSpace, dims_out, NULL);	int *Element_Connectivity; Element_Connectivity = (int*) malloc((dims_out[0]) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Connectivity); status=H5Sclose(DataSpace); status=H5Dclose(DataSet);

	// ///////////////////////////////////////////////// Reading Nodes dataset //////////////////////////////////////////////////////////////////////////////////////////////

	DataSet = H5Dopen(File, "Model/Nodes/Coordinates", H5P_DEFAULT); double *Node_Coordinates; Node_Coordinates = (double*) malloc((Number_of_Nodes*3) * sizeof(double));
	status = H5Dread(DataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,Node_Coordinates); status=H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Model/Nodes/Index_to_Coordinates", H5P_DEFAULT); int *Node_Index_to_Coordinates; Node_Index_to_Coordinates = (int*) malloc((Pseudo_Number_of_Nodes) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Node_Index_to_Coordinates); status=H5Dclose(DataSet);

	/////////////////////////////////////////////////// Building Nodes ///////////////////////////////////////////////////////////////////////////////////////////

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	points->SetNumberOfPoints(Number_of_Nodes);
	Node_Mesh->Allocate(Number_of_Elements);
	int node_no ;
	for (int i = 0; i < Number_of_Nodes; i++){
		node_no = Node_Map[i];
		points->InsertPoint(i, 
			Node_Coordinates[Node_Index_to_Coordinates[node_no]],
			Node_Coordinates[Node_Index_to_Coordinates[node_no]+1],
			Node_Coordinates[Node_Index_to_Coordinates[node_no]+2]
		);
	}

	Node_Mesh->SetPoints(points);

	///////////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	int connectivity_index,No_of_Element_Nodes,Cell_Type,element_no;

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

		Node_Mesh->InsertNextCell(Cell_Type, No_of_Element_Nodes, Vertices);

	}

	Whether_Node_Mesh_Build=1;

	status = H5Fclose(File);

	// freeing heap allocation for variables 
	free(Node_Map);
	free(Inverse_Node_Map);
	free(Element_Map);
	free(Element_Index_to_Connectivity);
	free(Element_Class_Tags);
	free(Element_Number_of_Nodes);
	free(Element_Connectivity);
	free(Node_Coordinates);
	free(Node_Index_to_Coordinates);

	return;
	
}

/*****************************************************************************
* This Function builds a gauss mesh from the given gauss mesh coordinates  
* Uses Delaunay3D filter to generate the mesh 
*****************************************************************************/

void pvESSI::Get_Gauss_Mesh(vtkSmartPointer<vtkUnstructuredGrid> Gauss_Mesh){

	herr_t status;
	hid_t File, DataSet, DataSpace; 
	hsize_t  dims_out[1];

	//////////////////////////////////////////////////////// Reading Hdf5 File ///////////////////////////////////////////////////////////////////////////////////////

	File = H5Fopen(this->FileName, H5F_ACC_RDWR, H5P_DEFAULT);

	//////////////////////////////////////////////////// Reading Map Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet = H5Dopen(File, "Maps/ElementMap", H5P_DEFAULT); int *Element_Map; Element_Map = (int*) malloc((Number_of_Elements) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Map); status=H5Dclose(DataSet);

	/////////////////////////////////////////////////////////////////////// Reading Element Data ///////////////////////////////////////////////////////////////////////////////////////////

	DataSet = H5Dopen(File, "Model/Elements/Index_to_Gauss_Point_Coordinates", H5P_DEFAULT); int *Element_Index_to_Gauss_Point_Coordinates; Element_Index_to_Gauss_Point_Coordinates = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Index_to_Gauss_Point_Coordinates); status=H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Model/Elements/Number_of_Gauss_Points", H5P_DEFAULT); int *Element_Number_of_Gauss_Points; Element_Number_of_Gauss_Points = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int));
	status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Number_of_Gauss_Points); status=H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Model/Elements/Gauss_Point_Coordinates", H5P_DEFAULT); DataSpace = H5Dget_space(DataSet);
	status  = H5Sget_simple_extent_dims(DataSpace, dims_out, NULL);	float *Element_Gauss_Point_Coordinates; Element_Gauss_Point_Coordinates = (float*) malloc((dims_out[0]) * sizeof(float));
	status = H5Dread(DataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Gauss_Point_Coordinates); status=H5Sclose(DataSpace); status=H5Dclose(DataSet);

	////////////////////////////////////////////////////////////////////// Building up the elements //////////////////////////////////////////////////////

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(Number_of_Gauss_Nodes); int Index_of_Gauss_Nodes=0, element_no;
	int No_of_Element_Gauss_Nodes; vtkIdType onevertex;

	for (int i=0; i < Number_of_Elements; i++){
		element_no = Element_Map[i];
		No_of_Element_Gauss_Nodes = Element_Number_of_Gauss_Points[element_no];
		for(int j=0; j<No_of_Element_Gauss_Nodes ; j++){
			points->InsertPoint(Index_of_Gauss_Nodes, Element_Gauss_Point_Coordinates[Element_Index_to_Gauss_Point_Coordinates[element_no]+3*j],Element_Gauss_Point_Coordinates[Element_Index_to_Gauss_Point_Coordinates[element_no]+3*j+1],Element_Gauss_Point_Coordinates[Element_Index_to_Gauss_Point_Coordinates[element_no]+3*j+2]);
			/// Adding 1D point to mesh 
			onevertex = Index_of_Gauss_Nodes;
			Gauss_Mesh->InsertNextCell(VTK_VERTEX, 1, &onevertex);

			/// updating index
			Index_of_Gauss_Nodes = Index_of_Gauss_Nodes+1;
		}

	}

	status = H5Fclose(File);

	Gauss_Mesh->SetPoints(points);

	this->Build_Delaunay3D_Gauss_Mesh(Gauss_Mesh);

	Whether_Gauss_Mesh_Build=1;

	// freeing heap allocation for variables  
	free(Element_Map);
	free(Element_Index_to_Gauss_Point_Coordinates);
	free(Element_Number_of_Gauss_Points);
	free(Element_Gauss_Point_Coordinates);

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
	}

	return;

}


/*****************************************************************************
* Initializes some variables to default values
*****************************************************************************/

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


/*****************************************************************************
* Builds maps in hdf5 output file for efficient visualization
* Creates a maps group in the files with the following data set 
*	Node_Map: From global/(input) node number to reduced node numbers
*	Element_Map : From global(input) element number to reduced element numbers
* 	Inverse_Node_Map: From reduced node number to global(input) node number 
* 	Inverse_Element_Map: From reduced element number to global(input) element number 
*****************************************************************************/

void pvESSI::Build_Maps(){

	herr_t status;
	hid_t File, DataSpace, DataSet, Group; 
	hsize_t  dims_out[1];

	//////////////////////////////////////////////////////// Reading Hdf5 File ///////////////////////////////////////////////////////////////////////////////////////

	File = H5Fopen(this->FileName, H5F_ACC_RDWR, H5P_DEFAULT);

	////////////////////////////////////////////////////// Reading General Outside Data /////////////////////////////////////////////////////////////////////////////

	DataSet = H5Dopen(File, "Model/Nodes/Index_to_Coordinates", H5P_DEFAULT); DataSpace = H5Dget_space(DataSet);
	status  = H5Sget_simple_extent_dims(DataSpace, dims_out, NULL);this->Pseudo_Number_of_Nodes = dims_out[0];
	int *Node_Index_to_Coordinates; Node_Index_to_Coordinates = (int*) malloc((Pseudo_Number_of_Nodes) * sizeof(int)); 
    status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT, Node_Index_to_Coordinates);
    H5Sclose(DataSpace); H5Dclose(DataSet);

	DataSet = H5Dopen(File, "Model/Elements/Index_to_Connectivity", H5P_DEFAULT); DataSpace = H5Dget_space(DataSet);
	status  = H5Sget_simple_extent_dims(DataSpace, dims_out, NULL);this->Pseudo_Number_of_Elements = dims_out[0];
	int *Element_Index_to_Connectivity; Element_Index_to_Connectivity = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int)); 
    status = H5Dread(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT, Element_Index_to_Connectivity);
    H5Sclose(DataSpace); H5Dclose(DataSet);

	int *Node_Map; Node_Map = (int*) malloc((Number_of_Nodes) * sizeof(int));
	int *Element_Map; Element_Map = (int*) malloc((Number_of_Elements) * sizeof(int));
	int *Inverse_Node_Map; Inverse_Node_Map = (int*) malloc((Pseudo_Number_of_Nodes) * sizeof(int));
	int *Inverse_Element_Map; Inverse_Element_Map = (int*) malloc((Pseudo_Number_of_Elements) * sizeof(int));
	int Whether_Maps_Build[1] = {1};

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

	//////////////////////// Debugg Printing /////////////////////////////////////////////////////
	// cout << "Number of Nodes " << index << " " << Number_of_Nodes << endl;
	// cout << "Pseudo_Number_of_Nodes " << index << " " << Pseudo_Number_of_Nodes << endl;
	// cout << "Number of Elements " << index << " " << Number_of_Elements << endl;;
	// cout << "Pseudo_Number_of_Elements " << index << " " << Pseudo_Number_of_Elements << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////

	Group = H5Gcreate(File, "/Maps", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); status = H5Gclose(Group); 

	dims_out[0]= Number_of_Nodes; DataSpace = H5Screate_simple(1, dims_out, NULL);
	DataSet = H5Dcreate(File,"Maps/NodeMap",H5T_STD_I32BE,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Node_Map);
	status = H5Sclose(DataSpace); status = H5Dclose(DataSet); 

	dims_out[0]= Number_of_Elements; DataSpace = H5Screate_simple(1, dims_out, NULL);
	DataSet = H5Dcreate(File,"Maps/ElementMap",H5T_STD_I32BE,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Element_Map);
	status = H5Sclose(DataSpace); status = H5Dclose(DataSet); 

	dims_out[0]= Pseudo_Number_of_Nodes; DataSpace = H5Screate_simple(1, dims_out, NULL);
	DataSet = H5Dcreate(File,"Maps/InverseNodeMap",H5T_STD_I32BE,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Inverse_Node_Map);
	status = H5Sclose(DataSpace); status = H5Dclose(DataSet); 

	dims_out[0]= Pseudo_Number_of_Elements; DataSpace = H5Screate_simple(1, dims_out, NULL);
	DataSet = H5Dcreate(File,"Maps/InverseElementMap",H5T_STD_I32BE,DataSpace,H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Inverse_Element_Map);
	status = H5Sclose(DataSpace); status = H5Dclose(DataSet); 

	DataSet = H5Dopen(File,"Whether_Maps_Build",H5P_DEFAULT);DataSpace = H5Dget_space(DataSet);
	status = H5Dwrite(DataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,Whether_Maps_Build);
	status = H5Sclose(DataSpace); status = H5Dclose(DataSet); 

	H5Fclose(File);

	// freeing heap allocation for variables 
	free(Node_Index_to_Coordinates);
	free(Element_Index_to_Connectivity);
}


/*****************************************************************************
* Builds time map i.e. a way to get the index number from a given non-interger
* time request from paraview desk. Although it works fine, its not the great 
* way to achieve it. 
* Needs to be changed. I think, paraview has fixed this issue in their latest 
* version, but I am not sure, need to check it. 
*****************************************************************************/

void pvESSI::Build_Time_Map(){

	for (int i = 0;i<Number_Of_Time_Steps ; i++)
		Time_Map[Time[i]] = i;

	return;
}


/*****************************************************************************
* Builds Meta Array map, so that it is easie to add attributes to vtkObject 
* and is neat and clean.
*****************************************************************************/

void pvESSI::Build_Meta_Array_Map(){

	int key = 0;

	Meta_Array_Map["Generalized_Displacements"] = key; key=key+1;
	Meta_Array_Map["Generalized_Velocity"] = key; key=key+1;
	Meta_Array_Map["Generalized_Acceleration"] = key; key=key+1;
	Meta_Array_Map["Elastic_Strain"] = key; key=key+1;
	Meta_Array_Map["Plastic_Strain"] = key; key=key+1;
	Meta_Array_Map["Stress"] = key; key=key+1;
	Meta_Array_Map["Material_Tag"] = key; key=key+1;
	Meta_Array_Map["Total_Energy"] = key; key=key+1;
	Meta_Array_Map["Incremental_Energy"] = key; key=key+1;
	Meta_Array_Map["Node_Tag"] = key; key=key+1;
	Meta_Array_Map["Element_Tag"] = key; key=key+1;

}



