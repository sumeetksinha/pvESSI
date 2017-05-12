pvESSI
=========

### Usage

It is plugin for visualizing the output of [REAL ESSI](http://sokocalo.engr.ucdavis.edu/~jeremic/Real_ESSI_Simulator/), a 3D finite element program specifically developed for high fidelity modeling and simulation of Earthquake Soil/Rock Structure Interaction problems. It reads ESSI (HDF5) .feioutput file.

### Installation

Currently, the plugin is not distributable type. There are two ways to install the plugin. First (A), is to build the plugin along with the paraview source and the other method (B) is to build the plugin externally and then integrate with paraview. The second option assumes that paraview is already build fro source into the system. For both ways, the steps are documented below. Although both ways require paraview to be present, it is encouraged to build pvESSI along with paraview source.

#### (A) Building Paraview from Source along with pvESSI

1. Download [ParaView source](http://www.paraview.org/download/) from its website and extract the folder.
2. Check or install these dependencies
	```bash
	sudo apt-get install libphonon-dev libphonon4 qt4-dev-tools 
	sudo apt-get install libqt4-core libqt4-gui qt4-qmake libxt-dev 
	sudo apt-get install g++ gcc cmake-curses-gui libqt4-opengl-dev 
	sudo apt-get install mesa-common-dev python-dev
	sudo apt-get install libavdevice-dev libavformat-dev libavfilter-dev libavcodec-dev libswscale-dev libavutil-dev
	sudo apt-get install ffmpeg # works for Ubuntu 16.04 or highger
	```
3. Check if you have the required cmake version for paraview source. Look in file **/Paraview/CMakeLists.txt** and also check your version. An example is shown for ParaView5.1.2

	```bash
	$ sed '31q;d' CMakeLists.txt
	cmake_minimum_required(VERSION 3.3)
	$ cmake --version
	cmake version 3.5.1
	```

  	If you have an older version, then get and install a newer version from [CMake](https://cmake.org/download/). You can also follow these steps at your desired directory.    
  	```bash
  	git clone https://github.com/Kitware/CMake.git
	cd CMake 
	mkdir build 
	cd build
	cmake ..
	make -j 8
	sudo make install
	```

4. Assuming, all the dependencies are satisfied, proceed with the installation in any directory  [${ParaView_Build_Directory}] outside ParaView source directory [${ParaView_Source_Directory}]. 

	(a) To use paraview in parallel build with ```PARAVIEW_USE_MPI=true``` option
	(b) If FFMPEG is already build in system with all it's libraries, build paraview with ```PARAVIEW_USE_MPI=true``` option
	```bash
	cd ${ParaView_Source_Directory}/Plugins/ 
	git clone https://github.com/SumeetSinha/pvESSI.git
	cd ${ParaView_Build_Directory}
	cmake -DPARAVIEW_USE_MPI=true -DPARAVIEW_ENABLE_PYTHON=true -DPARAVIEW_ENABLE_FFMPEG=true ${ParaView_Source_Directory}
	cmake ${ParaView_Source_Directory}
	make -j ${nop}
	```

	Where `$(nop)` is the number of processes you wish yo use.

5. If all goes well, ParaView gets build in the ${build_directory/bin}.
5. Now, you should be able to load and visualize ESSI (HDF5) .h5.feioutput file using paraview located at ```{ParaView_Build_Directory}/bin/paraview```.

#### (B) Compiling the pvESSI plugin externally into paraview

1. Get the latest version of plugin from github

	```bash
	git clone https://github.com/SumeetSinha/pvESSI.git
	```

2. Compile the plugin

	```bash
	cd pvESSI
	mkdir build
	cd ./build
	cmake .. -DParaView_DIR="${ParaView_Build_Directory}" 
	make -j ${nop}
	```

3. If the compilation suceeds, you will see ***libpvESSI.so*** built inside the build folder.
   Now, Open ParaView.

	```bash
	${ParaView_Build_Directory}/bin/paraview
	```

	a) Click on ***Tools*** and then on ***Manage Plugins***
	b) Click on ***Load New***. Navigate to the ***libpvESSI.so*** library and add it. 
	c) Check on the ***Auto Load*** option.

4. Close paraview and then agin open. Now, you should be able to load and visualize ESSI (HDF5) ``` .h5.feioutput``` file.


----
[UCD CompGeoMech](http://sokocalo.engr.ucdavis.edu/~jeremic/)

Created by: [Sumeet Kumar Sinha](http://www.sumeetksinha.com)

Request for features at sumeet.kumar507@gmail.com
   




