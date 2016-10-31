cmake_minimum_required(VERSION 2.8)

include_directories(/usr/lib/)

#add_subdirectory(Filter/Filter)

FIND_PACKAGE(ParaView REQUIRED)
LINK_LIBRARIES(hdf5 hdf5_cpp boost_regex  boost_system)
INCLUDE(${PARAVIEW_USE_FILE})
 ADD_PARAVIEW_PLUGIN(pvESSI "1.0"
   SERVER_MANAGER_XML pvESSI.xml
   SERVER_MANAGER_SOURCES pvESSI.cxx 
   GUI_RESOURCES pvESSI.qrc)
