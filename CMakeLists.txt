cmake_minimum_required (VERSION 2.8)
project (m)

# YCM config
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -fbuiltin -g")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
IF(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
                                                                ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json)
ENDIF()

file(GLOB sources "src/*.cpp")
list(REMOVE_ITEM sources "src/main.cpp")

add_library(m ${sources})
add_executable(testing "src/main.cpp")

include_directories(m "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(testing m)

install(TARGETS m DESTINATION "/usr/lib")

