cmake_minimum_required (VERSION 3.1)
project (mth)

option(TESTS "Whether to compile test program" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")

  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json)

endif()

# Build library

file(GLOB MTH_SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM MTH_SOURCES "${PROJECT_SOURCE_DIR}/src/main.cpp")

add_library(mth STATIC ${MTH_SOURCES})
target_include_directories(mth PUBLIC "${PROJECT_SOURCE_DIR}/include")

# Build tests

if(TESTS)

    configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)

    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .

        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download
    )

    if(result)

        message(FATAL_ERROR "CMake step for googletest failed: ${result}")

    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} --build .

        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download
    )

    if(result)

        message(FATAL_ERROR "Build step for googletest failed: ${result}")

    endif()

    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src

        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL
    )

    if(CMAKE_VERSION VERSION_LESS 2.8.11)
        
        include_directories("${gtest_SOURCE_DIR}/include")

    endif()

    add_executable(mth_test "${PROJECT_SOURCE_DIR}/src/main.cpp")
    target_link_libraries(mth_test mth gtest_main)

    add_test(NAME tests COMMAND mth_test)
endif()
