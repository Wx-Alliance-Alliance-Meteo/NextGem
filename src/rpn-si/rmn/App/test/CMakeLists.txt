if (EC_INIT_DONE LESS 2)
   include(CTest)
   add_custom_target(check COMMAND CTEST_OUTPUT_ON_FAILURE=true ${CMAKE_CTEST_COMMAND})

   if(BUILD_TESTING)
        set(CMAKE_Fortran_PREPROCESS ON)
        add_executable(test_app_loglevel_noinit EXCLUDE_FROM_ALL test_app_loglevel_noinit.f90)
        add_test(
            NAME test_app_loglevel_noinit
            COMMAND $<TARGET_FILE:test_app_loglevel_noinit>
        )
        set_tests_properties(test_app_loglevel_noinit
            PROPERTIES PASS_REGULAR_EXPRESSION "INFO")
        target_include_directories(test_app_loglevel_noinit PRIVATE ${CMAKE_BINARY_DIR}/src/include)
        target_link_libraries(test_app_loglevel_noinit App::App)
        add_dependencies(check test_app_loglevel_noinit)
    endif()
endif()
