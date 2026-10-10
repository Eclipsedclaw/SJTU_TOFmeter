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
# Qt helpers for the GUI when Geant4 was built with Qt: /tof/vis/spin, and the TOFmeter logo as
# the application icon (the marks of the repository's data/ folder, embedded in the executable;
# needs the Qt SVG module)
if(Geant4_qt_FOUND)
  if(TARGET Qt6::Widgets)
    set(TOF_QT Qt6)
  elseif(TARGET Qt5::Widgets)
    set(TOF_QT Qt5)
  endif()
  if(TOF_QT)
    target_link_libraries(TOFmeter ${TOF_QT}::Widgets)
    target_compile_definitions(TOFmeter PRIVATE TOF_USE_QT)
    find_package(${TOF_QT} QUIET COMPONENTS Svg)
    set(TOF_LOGO_MARK_FILE ${PROJECT_SOURCE_DIR}/../data/tofmeter-logo-mark.svg)
    set(TOF_LOGO_REVERSED_FILE ${PROJECT_SOURCE_DIR}/../data/tofmeter-logo-mark-reversed.svg)
    if(TARGET ${TOF_QT}::Svg AND EXISTS ${TOF_LOGO_MARK_FILE} AND EXISTS ${TOF_LOGO_REVERSED_FILE})
      file(READ ${TOF_LOGO_MARK_FILE} TOF_LOGO_MARK)
      file(READ ${TOF_LOGO_REVERSED_FILE} TOF_LOGO_MARK_REVERSED)
      configure_file(${PROJECT_SOURCE_DIR}/include/TOFmeterLogo.hh.in ${PROJECT_BINARY_DIR}/TOFmeterLogo.hh @ONLY)
      # re-run CMake when a logo changes
      set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${TOF_LOGO_MARK_FILE} ${TOF_LOGO_REVERSED_FILE})
      target_include_directories(TOFmeter PRIVATE ${PROJECT_BINARY_DIR})
      target_link_libraries(TOFmeter ${TOF_QT}::Svg)
      target_compile_definitions(TOFmeter PRIVATE TOF_HAS_LOGO)
    else()
      message(STATUS "TOFmeter: no logo icon (needs ${TOF_QT}Svg and data/tofmeter-logo-mark*.svg)")
    endif()
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
