if (EC_INIT_DONE LESS 2)
   if(BUILD_TESTING)

      #----- Library tests
      message(STATUS "(EC) Generating vgrid library and utilities tests")

      file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/Testing/data_out ${CMAKE_BINARY_DIR}/Testing/data)
      execute_process(COMMAND ln -s $ENV{ECCI_DATA_DIR}/vgrid/b_6.5/data ${CMAKE_BINARY_DIR}/Testing/data)
      execute_process(COMMAND rm -f -r ${CMAKE_BINARY_DIR}/Testing/data_out)
      execute_process(COMMAND mkdir ${CMAKE_BINARY_DIR}/Testing/data_out)

      file(GLOB SRC_TEST lib/src/*.c lib/src/*.F90)

      #----- Loop on tests
      foreach(testfile ${SRC_TEST})
         get_filename_component(ext  ${testfile} EXT)
         get_filename_component(test ${testfile} NAME)
         get_filename_component(test ${test}     NAME_WLE)
   #      message(STATUS "(EC) Adding ${ext} test ${testfile}")

         if(${ext} MATCHES ".c")
            add_executable(${test} EXCLUDE_FROM_ALL ${testfile} lib/c_ut_report.c)
         else()
            add_executable(${test} EXCLUDE_FROM_ALL ${testfile} lib/unit_testing.F90)
         endif()
         add_dependencies(check ${test})

         target_link_libraries(${test} PRIVATE vgrid::vgrid)
         target_include_directories(${test} PRIVATE ${CMAKE_BINARY_DIR}/lib ${CMAKE_SOURCE_DIR}/src/lib ${CMAKE_CURRENT_SOURCE_DIR} ${RMN_INCLUDE_DIRS})

         add_test(NAME ${test} COMMAND ${test} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing)
      endforeach()

      #----- Utilities tests
      find_program(BASH_COMMAND bash)
      add_test(NAME testutils COMMAND ${BASH_COMMAND} ${CMAKE_SOURCE_DIR}/test/utils/testutilvgrid.sh ${CMAKE_BINARY_DIR}/utils $ENV{ECCI_DATA_DIR}/vgrid/1.9.b3 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing)
   endif()
endif()
