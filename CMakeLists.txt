# CMakeLists.txt
cmake_minimum_required(VERSION 3.0)
project(framework CXX)

# Find dependencies provided by vcpkg (via vcpkg.cmake)
find_package(glfw3 CONFIG REQUIRED)
find_package(glm CONFIG REQUIRED)
find_package(assimp CONFIG REQUIRED)
find_package(glad CONFIG REQUIRED)
find_package (Eigen3 3.3 NO_MODULE)
find_package(spdlog REQUIRED)
set(PROJECT_NAME main)

################################################################################
# Source groups
################################################################################
set(Headers
    "src/camera.h"
    "src/light.h"
    "src/mesh.h"
    "src/model.h"
    "src/shader.h"
    "src/stb_image.h"
    "src/utils.h"
    "src/dynamics.h"
    "src/globals.h"
)
source_group("Headers" FILES ${Headers})

set(Sources
    "src/light.cpp"
    "src/main.cpp"
    "src/dynamics.cpp"
)
source_group("Sources" FILES ${Sources})

set(ALL_FILES
    ${Headers}
    ${Sources}
)

################################################################################
# Target
################################################################################
add_executable(${PROJECT_NAME} ${ALL_FILES})

# # main target
# add_executable(main)
target_include_directories(main PRIVATE 
    ${GLM_INCLUDE_DIRS} 
    ${GLFW3_INCLUDE_DIRS} 
    ${ASSIMP_INCLUDE_DIRS} 
    ${GLAD_INCLUDE_DIRS})

target_sources(main PRIVATE ${ALL_FILES})
target_link_libraries(main PRIVATE glfw assimp::assimp glad::glad glm::glm Eigen3::Eigen spdlog::spdlog_header_only)
# set_property(TARGET main PROPERTY CXX_STANDARD 20)

# # tests target
# add_executable(tests)
# target_sources(tests PRIVATE src/main.cpp)
# target_link_libraries(tests PRIVATE unofficial::sqlite3::sqlite3 fmt::fmt)
# set_property(TARGET tests PROPERTY CXX_STANDARD 20)
# enable_testing()
# add_test(tests tests)
