message(STATUS "Generating App library")

file(GLOB PROJECT_INCLUDE_FILES *.h)
file(GLOB PROJECT_C_FILES *.c)
file(GLOB PROJECT_F_FILES *.f90)

#----- Non ompi version
set(targets App App-shared)
set(shared_targets App-shared)
set(static_targets App-static)

add_library(App-static STATIC ${PROJECT_INCLUDE_FILES} ${PROJECT_C_FILES} ${PROJECT_F_FILES})
target_include_directories(App-static PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_include_directories(App-static INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>)

add_library(App-shared SHARED $<TARGET_OBJECTS:App-static>)
target_include_directories(App-shared PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_include_directories(App-shared INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>)
 
add_dependencies(App-static ${PROJECT_NAME}_build_info)
set_target_properties(App-static App-shared PROPERTIES
   VERSION ${PROJECT_VERSION}
   PUBLIC_HEADER "${PROJECT_INCLUDE_FILES}"
   INTERFACE_COMPILE_DEFINITIONS HAVE_APP="${PROJECT_VERSION}"
   Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include
   OUTPUT_NAME App
   POSITION_INDEPENDENT_CODE  ON)

#----- Needed when used as sbumodule/subtree
add_library(App::App-shared ALIAS App-shared)
add_library(App::App-static ALIAS App-static)
if(BUILD_SHARED_LIBS)
    add_library(App::App ALIAS App-shared)
else()
    add_library(App::App ALIAS App-static)
endif()

if(WITH_OMPI AND OpenMP_FOUND AND MPI_FOUND) 
   message(STATUS "Generating App ompi library")

   #----- ompi version
   list(APPEND targets App-ompi-static App-ompi-shared)
   list(APPEND shared_targets App-ompi-shared)
   list(APPEND static_targets App-ompi-static)

   add_library(App-ompi-static STATIC ${PROJECT_INCLUDE_FILES} ${PROJECT_C_FILES} ${PROJECT_F_FILES})
   add_dependencies(App-ompi-static ${PROJECT_NAME}_build_info)
   target_include_directories(App-ompi-static PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
   target_include_directories(App-ompi-static INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/ompi>)
   target_compile_definitions(App-ompi-static PUBLIC HAVE_MPI HAVE_OPENMP)
   target_link_libraries(App-ompi-static PUBLIC MPI::MPI_C OpenMP::OpenMP_C)

   add_library(App-ompi-shared SHARED $<TARGET_OBJECTS:App-ompi-static>)
   target_include_directories(App-ompi-shared PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
   target_include_directories(App-ompi-shared INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/ompi>)
   target_compile_definitions(App-ompi-shared PUBLIC HAVE_MPI HAVE_OPENMP)
   target_link_libraries(App-ompi-shared PUBLIC MPI::MPI_C OpenMP::OpenMP_C)

   set_target_properties(App-ompi-static App-ompi-shared PROPERTIES
      VERSION ${PROJECT_VERSION}
      PUBLIC_HEADER "${PROJECT_INCLUDE_FILES}"
      INTERFACE_COMPILE_DEFINITIONS HAVE_APP="${PROJECT_VERSION}"
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/ompi
      OUTPUT_NAME App-ompi
      POSITION_INDEPENDENT_CODE  ON)

   #----- Needed when used as sbumodule/subtree
   add_library(App::App-ompi-shared ALIAS App-ompi-shared)
   add_library(App::App-ompi-static ALIAS App-ompi-static)
      if(BUILD_SHARED_LIBS)
      add_library(App::App-ompi ALIAS App-shared)
   else()
      add_library(App::App-ompi ALIAS App-static)
   endif()
endif()

install(TARGETS ${shared_targets}
   EXPORT ${PROJECT_NAME}-shared-targets
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
   INCLUDES DESTINATION include
   PUBLIC_HEADER DESTINATION include)

install(TARGETS ${static_targets}
   EXPORT ${PROJECT_NAME}-static-targets
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib
   INCLUDES DESTINATION include
   PUBLIC_HEADER DESTINATION include)  

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/app.mod TYPE INCLUDE)
