cmake_minimum_required(VERSION 3.10)

# set the project name
project(projet-mnt)

# add the executable
add_executable(output main.cpp)

find_package(PROJ)

target_link_libraries(output PRIVATE ${PROJ_LIBRARIES})
