cmake_minimum_required(VERSION 3.10)

set(CMAKE_BUILD_TYPE Release)

# set the project name
project(projet-mnt)

# add the executable
add_executable(output main.cpp)


find_package(PROJ)
find_package(CGAL)
target_link_libraries(output ${PROJ_LIBRARIES} CGAL::CGAL)
