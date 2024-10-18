cmake_minimum_required(VERSION 3.20..3.25)
project(EVP 
    VERSION 1.0.0
    DESCRIPTION "Eigen value problem."
    )
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Testing
add_subdirectory(${PROJECT_SOURCE_DIR}/external)
include(CTest)
add_subdirectory(${PROJECT_SOURCE_DIR}/tests)

set(SRC_FILES ${PROJECT_SOURCE_DIR}/scripts/main.cxx)
add_executable(${PROJECT_NAME} ${SRC_FILES})
add_subdirectory(${PROJECT_SOURCE_DIR}/src)
list(APPEND INC_LIBS
    Linalg
)
target_link_libraries(${PROJECT_NAME} PRIVATE ${INC_LIBS})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)
