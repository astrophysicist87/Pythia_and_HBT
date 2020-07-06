# Install script for directory: /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/fio

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/HepMC/examples/fio" TYPE FILE FILES
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/fio/example_MyHerwig.cc"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/fio/example_MyPythia.cc"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/fio/example_MyPythiaOnlyToHepMC.cc"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/fio/example_PythiaStreamIO.cc"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/fio/initpydata.f"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/fio/initPythia.cc"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/fio/PythiaHelper.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/fio/testHerwigCopies.cc"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/fio/testPythiaCopies.cc"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/examples/fio/GNUmakefile"
    )
endif()

