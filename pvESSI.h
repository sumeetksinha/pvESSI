#ifndef __pvESSI_h
#define __pvESSI_h
 
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkUnstructuredGridAlgorithm.h"
#include <map>
#include <vector>
 
class pvESSI : public vtkUnstructuredGridAlgorithm
{
public:
  vtkTypeMacro(pvESSI,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static pvESSI *New();
 
  // Description:
  // Specify file name of the .abc file.
  vtkSetStringMacro(FileName);
  // vtkGetStringMacro(FileName);
 
protected:
  pvESSI();
  ~pvESSI(){}
 
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);
 
private:
  std::map<int,int> ESSI_to_VTK_Element;
  std::map< int,std::vector<int> > ESSI_to_VTK_Connectivity;
  pvESSI(const pvESSI&);  // Not implemented.
  void operator=(const pvESSI&);  // Not implemented.
  void set_VTK_To_ESSI_Elements_Connectivity();
 
  char* FileName;
};
 
#endif
