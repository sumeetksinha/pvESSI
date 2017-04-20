cmake_minimum_required(VERSION 2.8)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 4.9
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        message(FATAL_ERROR "GCC version must be at least 4.9!")
    endif() 
else()
    message(WARNING "You are using an unsupported compiler! Compilation has only been tested with GCC.")
endif()

IF (ParaView_SOURCE_DIR)
  INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIRS}
  )
ELSE()
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF ()

#include_directories(/usr/lib/)
#LINK_LIBRARIES(hdf5 hdf5_cpp boost_regex  boost_system)

include_directories(${VTK_INCLUDE_DIRS})
#include(wrap-python.cmake)

ADD_PARAVIEW_PLUGIN(pvESSI "1.0"
   SERVER_MANAGER_XML pvESSI.xml
   SERVER_MANAGER_SOURCES pvESSI.cxx)

target_link_libraries(pvESSI ${deps})

#if(VTK_WRAP_PYTHON)
#       wrap_python(pvESSI "${sources}")
#endif()

#vtk_module_library(vtkPvESSI ${Module_SRCS})


#from paraview.simple import *
#LoadPlugin("libpvESSI.so",ns=globals())
#from paraview.simple import *
#pvESSI(FileName="/home/sumeet/Desktop/ShearBox_Sequential.feioutput")
#Show()
#Render()
#viewLayout = GetLayout()
                          
