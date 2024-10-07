cmake_minimum_required(VERSION 3.20)
project(EVP VERSION 1.0.0)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(SRC_FILES ${PROJECT_SOURCE_DIR}/scripts/main.cxx)
add_executable(${PROJECT_NAME} ${SRC_FILES})
add_subdirectory(${PROJECT_SOURCE_DIR}/src)
list(APPEND INC_LIBS
    Linalg
)
target_link_libraries(${PROJECT_NAME} PRIVATE ${INC_LIBS})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)
