cmake_minimum_required(VERSION 3.20..3.25)
project(EVP 
    VERSION 1.0.0
    DESCRIPTION "Eigen value problem in parallel."
    )
option(WITH_MPI "Use MPI" ON)
if(WITH_MPI)
    find_package(MPI COMPONENTS CXX)
    if(MPI_FOUND)
        set(CMAKE_CXX_COMPILER "/usr/bin/mpicxx")
        set(SRC_FILES ${PROJECT_SOURCE_DIR}/scripts/mpi_script.cxx)
    else()
        set(SRC_FILES ${PROJECT_SOURCE_DIR}/scripts/sequential.cxx)
        message("-- MPI not found, using default C++ compiler")
        option(WITH_MPI OFF)
    endif(MPI_FOUND)
else()
    set(SRC_FILES ${PROJECT_SOURCE_DIR}/scripts/sequential.cxx)
    message("-- MPI disabled, using default C++ compiler")
endif(WITH_MPI)
# Update the submodules here
include(cmake/updateSubmodules.cmake)
# Testing
option(BUILD_TESTING ~WITH_MPI)
if(BUILD_TESTING)
    include(CTest)
    add_subdirectory(${PROJECT_SOURCE_DIR}/external)
    add_subdirectory(${PROJECT_SOURCE_DIR}/tests)
endif(BUILD_TESTING)
# Source
add_executable(${PROJECT_NAME} ${SRC_FILES})
add_subdirectory(${PROJECT_SOURCE_DIR}/src)
list(APPEND INC_LIBS
    Linalg
)
target_link_libraries(${PROJECT_NAME} PRIVATE ${INC_LIBS})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)
