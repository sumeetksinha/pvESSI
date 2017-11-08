cmake_minimum_required(VERSION 2.8)

#INCLUDE_DIRECTORIES(/home/sumeet/Softwares/RealESSI_ROOT/RealESSI_Dependencies/hdf5_swmr/include)
#LINK_DIRECTORIES(/home/sumeet/Softwares/RealESSI_ROOT/RealESSI_Dependencies/hdf5_swmr/lib)

# INCLUDE_DIRECTORIES(/home/sumeet/Softwares/RealESSI_ROOT/RealESSI_Dependencies/hdf5_sequential/include)
# LINK_DIRECTORIES(/home/sumeet/Softwares/RealESSI_ROOT/RealESSI_Dependencies/hdf5_sequential/lib)

# LINK_LIBRARIES(hdf5)
#LINK_LIBRARIES(hdf5_cpp hdf5 )

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




# include_directories(${VTK_INCLUDE_DIRS})
#include(wrap-python.cmake)

ADD_PARAVIEW_PLUGIN(pvESSI "1.0"
   SERVER_MANAGER_XML pvESSI.xml
   SERVER_MANAGER_SOURCES pvESSI.cxx)

# target_link_libraries(pvESSI hdf5 ${deps})
