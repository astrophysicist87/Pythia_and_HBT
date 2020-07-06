# Install script for directory: /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/HepMC" TYPE FILE FILES
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/CompareGenEvent.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/Flow.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/GenEvent.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/GenParticle.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/GenVertex.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/GenCrossSection.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/GenRanges.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/HeavyIon.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/HEPEVT_Wrapper.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/HerwigWrapper.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/IO_AsciiParticles.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/IO_BaseClass.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/IO_Exception.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/IO_GenEvent.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/IO_HEPEVT.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/IO_HERWIG.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/IteratorRange.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/PdfInfo.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/Polarization.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/PythiaWrapper6_4.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/PythiaWrapper6_4_WIN32.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/PythiaWrapper.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/WeightContainer.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/SearchVector.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/SimpleVector.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/SimpleVector.icc"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/StreamHelpers.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/StreamInfo.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/enable_if.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/is_arithmetic.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/TempParticleMap.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/Units.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/Version.h"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/HepMC/HepMCDefs.h"
    )
endif()

