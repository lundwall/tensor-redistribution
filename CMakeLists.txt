cmake_minimum_required(VERSION 3.9)
project(dphpc)
find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_BUILD_TYPE "Release")

add_executable(2dtest src/2d_distribute.cpp)
target_link_libraries(2dtest PUBLIC MPI::MPI_CXX)