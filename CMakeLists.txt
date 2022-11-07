cmake_minimum_required(VERSION 3.9)
project(dphpc)
find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(datatype)
add_subdirectory(manual)
add_subdirectory(codegen_sample)