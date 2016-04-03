cmake_minimum_required(VERSION 2.8)

include_directories(/usr/lib /usr/include)
FIND_PACKAGE(ParaView REQUIRED)
LINK_LIBRARIES(hdf5_cpp hdf5 )
INCLUDE(${PARAVIEW_USE_FILE})
 ADD_PARAVIEW_PLUGIN(pvESSI "1.0"
   SERVER_MANAGER_XML pvESSI.xml
   SERVER_MANAGER_SOURCES pvESSI.cxx 
   GUI_RESOURCES pvESSI.qrc)