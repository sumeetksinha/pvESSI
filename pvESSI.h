#ifndef __pvESSI_h
#define __pvESSI_h
 
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkUnstructuredGridAlgorithm.h"
#include <map>
#include <vector>
#include <iostream>

class pvESSI : public vtkUnstructuredGridAlgorithm
{
public:
  vtkTypeMacro(pvESSI,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  void PrintX(int x){
  	std::cout << x << std::endl;
  }
  // void GetNumberOfCellArrays();
 
  static pvESSI *New();
 
  // Description:
  // Specify file name of the .abc file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
  vtkSetMacro(TimeStep, int);
  vtkGetMacro(TimeStep, int);
  int GetNumberOfTimeSteps(){
    return this->NumberOfTimeSteps;
  }
 
protected:
  pvESSI();
  ~pvESSI(){}
 
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestInformation( vtkInformation *, vtkInformationVector **, vtkInformationVector* );
  static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);
  
  //////////////////////// Important Variables /////////////////////////////////////////

  int TimeStep;
  int CurrentTimeStep; 
  int NumberOfTimeSteps;
  vtkTimeStamp ReadMTime;
  int ReadError;
  int TimeStepRange[2];

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
  std::map< int,std::vector<int> > ESSI_to_VTK_Connectivity;
  pvESSI(const pvESSI&);  // Not implemented.
  void operator=(const pvESSI&);  // Not implemented.
  void set_VTK_To_ESSI_Elements_Connectivity();
 
  char* FileName;
};
 
#endif
