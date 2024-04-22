message(STATUS "(EC) Generating vgrid examples")

file(GLOB SRC_EXAMPLES *.c *.F90)

#----- Loop on tests
foreach(exfile ${SRC_EXAMPLES})
   get_filename_component(ex ${exfile} NAME)
   get_filename_component(ex ${exfile} NAME_WLE)

   add_executable(${ex} ${exfile})
   target_link_libraries(${ex} PRIVATE vgrid)
   target_include_directories(${ex} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}  ${CMAKE_SOURCE_DIR}/src/lib  ${rmn_INCLUDE_DIRS})
endforeach()
