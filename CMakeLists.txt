cmake_minimum_required (VERSION 2.8.2)
project (mth)

# YCM config
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -fbuiltin -g")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")

  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json)

endif()

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

file(GLOB sources "${CMAKE_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM sources "${CMAKE_SOURCE_DIR}/src/main.cpp")

add_library(mth_lib ${sources})
add_executable(mth_test "${CMAKE_SOURCE_DIR}/src/main.cpp")
add_test(NAME tests COMMAND mth_test)

include_directories(mth_lib "${CMAKE_SOURCE_DIR}/include")
include_directories(mth_test "${CMAKE_SOURCE_DIR}/include")

target_link_libraries(mth_test mth_lib gtest_main)

