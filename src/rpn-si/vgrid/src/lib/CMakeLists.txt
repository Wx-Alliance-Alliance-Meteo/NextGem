message(STATUS "(EC) Generating vgrid librairie")

set(PROJECT_INCLUDE_FILES vgrid.h)
set(MODULE_DIR ${CMAKE_CURRENT_BINARY_DIR}/include)
set(PROJECT_MODULE_FILES ${MODULE_DIR}/vgrid_descriptors.mod ${MODULE_DIR}/vgrid_utils.mod)

add_library(vgrid-static STATIC vgrid_descriptors.F90 vgrid_utils.F90 vgrid.c)
add_dependencies(vgrid-static ${PROJECT_NAME}_build_info)
add_library(vgrid-shared SHARED $<TARGET_OBJECTS:vgrid-static>)
add_dependencies(vgrid-shared vgrid-static ${PROJECT_NAME}_build_info)

set_target_properties(vgrid-static vgrid-shared PROPERTIES
  VERSION ${PROJECT_VERSION}
  PUBLIC_HEADER "${PROJECT_INCLUDE_FILES}"
  Fortran_MODULE_DIRECTORY ${MODULE_DIR}
  INTERFACE_COMPILE_DEFINITIONS HAVE_VGRID="${PROJECT_VERSION}"
  INTERFACE_INCLUDE_DIRECTORIES $<BUILD_INTERFACE:${MODULE_DIR}>
  POSITION_INDEPENDENT_CODE ON
  OUTPUT_NAME vgrid)
 
target_link_libraries(vgrid-static rmn::rmn-static)
target_link_libraries(vgrid-shared rmn::rmn-shared)

install(TARGETS vgrid-static
  EXPORT vgrid-static-targets
  INCLUDES      DESTINATION include
  PUBLIC_HEADER DESTINATION include
  LIBRARY       DESTINATION lib
  ARCHIVE       DESTINATION lib)
install(TARGETS vgrid-shared
  EXPORT vgrid-shared-targets
  INCLUDES      DESTINATION include
  PUBLIC_HEADER DESTINATION include
  LIBRARY       DESTINATION lib
  ARCHIVE       DESTINATION lib)
install(FILES ${PROJECT_MODULE_FILES} DESTINATION include)

if(BUILD_SHARED_LIBS)
    add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME}-shared)
else()
    add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME}-static)
endif()

# Symbolic link for archive numbered version
add_custom_command(TARGET ${PROJECT_NAME}-static POST_BUILD
    COMMAND ln -sfv lib${PROJECT_NAME}.a lib${PROJECT_NAME}_${PROJECT_VERSION}.a)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}_${PROJECT_VERSION}.a DESTINATION lib)
