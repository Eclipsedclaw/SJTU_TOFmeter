# @author Tsuguo Aramaki
# @date 2015 March 23

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(TOFmeter)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Optional: Geant4 dataset directory used when GEANT4_DATA_DIR is not set at
# run time, e.g. -DTOF_GEANT4_DATA_DIR=/path/to/geant4-build/data
#
set(TOF_GEANT4_DATA_DIR "" CACHE PATH "Fallback for the GEANT4_DATA_DIR environment variable")

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh
                  ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(TOFmeter TOFmeter.cc ${sources} ${headers})
target_compile_features(TOFmeter PRIVATE cxx_std_17)
target_compile_definitions(TOFmeter PRIVATE TOF_GEANT4_DATA_DIR="${TOF_GEANT4_DATA_DIR}")
target_link_libraries(TOFmeter ${Geant4_LIBRARIES})
# Qt helpers for the GUI (/tof/vis/spin) when Geant4 was built with Qt
if(Geant4_qt_FOUND)
  if(TARGET Qt6::Widgets)
    target_link_libraries(TOFmeter Qt6::Widgets)
    target_compile_definitions(TOFmeter PRIVATE TOF_USE_QT)
  elseif(TARGET Qt5::Widgets)
    target_link_libraries(TOFmeter Qt5::Widgets)
    target_compile_definitions(TOFmeter PRIVATE TOF_USE_QT)
  endif()
endif()

#----------------------------------------------------------------------------
# Copy all macros, analysis scripts and data tables to the build directory, i.e. the directory in which we
# build TOFmeter. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
file(GLOB TOFmeter_SCRIPTS RELATIVE ${PROJECT_SOURCE_DIR}
     ${PROJECT_SOURCE_DIR}/macro/*.mac
     ${PROJECT_SOURCE_DIR}/analysis/*.C
     ${PROJECT_SOURCE_DIR}/scripts/*.sh
     ${PROJECT_SOURCE_DIR}/data/*.txt)

foreach(_script ${TOFmeter_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(PJ DEPENDS TOFmeter)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS TOFmeter DESTINATION bin)
