cmake_minimum_required(VERSION 3.10)

set(CMAKE_BUILD_TYPE Release)

# set the project name
project(projet-mnt)

# add the executable
add_executable(create_raster main.cpp)

set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(PROJ)
find_package(CGAL)
find_package(Threads REQUIRED)
target_link_libraries(create_raster ${PROJ_LIBRARIES} CGAL::CGAL Threads::Threads)
