# Install script for directory: /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/pythia8

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/HepMC/examples/pythia8" TYPE FILE FILES
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/pythia8/main31.cc"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/pythia8/main32.cc"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/pythia8/main32.cmnd"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/examples/pythia8/README"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/examples/pythia8/GNUmakefile"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/examples/pythia8/config.csh"
    "/scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/build/examples/pythia8/config.sh"
    )
endif()

