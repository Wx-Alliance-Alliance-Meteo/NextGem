cmake_minimum_required(VERSION 3.16)

message(STATUS "(EC) Generating rpncomm Makefile")

#----- Append EC specific module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_rpn/modules $ENV{EC_CMAKE_MODULE_PATH})

include(ec_init)
#ec_git_version()           # Get version from git state
set(VERSION "4.6.2-8-g621992f")
ec_parse_manifest()

project(${NAME} VERSION ${VERSION} DESCRIPTION "${DESCRIPTION}") 
set(PROJECT_VERSION ${VERSION}${STATE})
set(rpncomm_VERSION ${PROJECT_VERSION} CACHE INTERNAL "rpncomm version" FORCE) # Needed for cascaded version identification
message(STATUS "(EC) ${PROJECT_NAME} version = ${PROJECT_VERSION}")

string(REPLACE "." "" RPN_COMM_VERSION_MARKER ${VERSION})

option(RPN_SERIAL "Whether the library will avoid MPI (not fully implemented, DO NOT USE)")

if (RPN_SERIAL)
   message(WARNING "(EC) The RPN_SERIAL option is not fully implemented. It will be deactivated.")
   set(RPN_SERIAL FALSE)
endif()

enable_language(Fortran)
enable_language(C)
include(ec_compiler_presets)

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

# Find out MPI version (to use in rpncomm version)
if (NOT RPN_SERIAL)
   find_program(RPN_MPIRUN mpirun)
   if (${RPN_MPIRUN} STREQUAL "RPN_MPIRUN-NOTFOUND")
      message(FATAL_ERROR "(EC) mpirun not found. Won't be able to configure the version.")
   else()
     execute_process(
       COMMAND bash "-c" "${RPN_MPIRUN} --version 2>&1 | head -1 | xargs -n1 | grep '\\.'"
       OUTPUT_VARIABLE MPI_LIBRARY_VERSION
       OUTPUT_STRIP_TRAILING_WHITESPACE)
     set(MPI_LIBRARY_VERSION "_${MPI_LIBRARY_VERSION}")
   endif()
endif()

add_subdirectory(src)

#----- Generate the config file for the project to be usable via cmake's find_package command
set(INCLUDE_INSTALL_DIR include)
set(LIB_INSTALL_DIR     lib)
set(CONFIG_INSTALL_DIR  "${LIB_INSTALL_DIR}/cmake/${PROJECT_NAME}-${PROJECT_VERSION}")

include(CMakePackageConfigHelpers)
configure_package_config_file(
    "Config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${CONFIG_INSTALL_DIR}"
    PATH_VARS           INCLUDE_INSTALL_DIR LIB_INSTALL_DIR
)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion
)

install(EXPORT ${PROJECT_NAME}-targets
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CONFIG_INSTALL_DIR})

install(FILES   "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION "${CONFIG_INSTALL_DIR}"
)

#----- Packaging
ec_package_name()    # Define package prefix  
ec_build_config()    # Create build configuration script
ec_prepare_ssm()     # Prepare ssm packaging files

set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "ECCC")
set(CPACK_PACKAGE_CONTACT "${MAINTAINER}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.org")
set(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/package")
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NAME}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${NAME}_${PROJECT_VERSION}")
include(CPack)

if("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${PACKAGE_NAME}" CACHE FILEPATH "CMake Installation prefix" FORCE)
endif()
