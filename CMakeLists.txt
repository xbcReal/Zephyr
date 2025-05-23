cmake_minimum_required(VERSION 3.10)

project(HAWK VERSION 1.0)

# CMAKE config
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

# Debugging
set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add dependencies
find_package(SEAL 4.1 REQUIRED)

# Nexus
add_subdirectory(src/bootstrapping/common)
add_subdirectory(src/bootstrapping)

add_executable(
    main
    ${CMAKE_SOURCE_DIR}/src/main.cpp
    ${CMAKE_SOURCE_DIR}/src/gelu.cpp
    ${CMAKE_SOURCE_DIR}/src/layer_norm.cpp
    ${CMAKE_SOURCE_DIR}/src/ckks_evaluator.cpp
    ${CMAKE_SOURCE_DIR}/src/softmax.cpp
    ${CMAKE_SOURCE_DIR}/src/matrix_mul.cpp
    ${CMAKE_SOURCE_DIR}/src/argmax.cpp
    ${COMMON_SOURCE_FILES}
    ${BOOTSTRAPPING_SOURCE_FILES}
)

target_include_directories(main PUBLIC
    /usr/local/include
    /usr/local/include/NTL
    ${COMMON_HEADER_DIR}
    ${BOOTSTRAPPING_HEADER_DIR}
)

target_link_libraries(main PRIVATE ntl gmp m pthread SEAL::seal)

# Nexus Bootstrapping
add_executable(
    bootstrapping
    ${CMAKE_SOURCE_DIR}/src/bootstrapping.cpp
    ${CMAKE_SOURCE_DIR}/src/ckks_evaluator.cpp
    ${COMMON_SOURCE_FILES}
    ${BOOTSTRAPPING_SOURCE_FILES}
)

target_include_directories(bootstrapping PUBLIC
    /usr/local/include
    /usr/local/include/NTL
    ${COMMON_HEADER_DIR}
    ${BOOTSTRAPPING_HEADER_DIR}
)

target_link_libraries(bootstrapping PRIVATE ntl gmp m pthread SEAL::seal)

