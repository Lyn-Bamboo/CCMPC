cmake_minimum_required(VERSION 2.8.3)
project(grid_map)

find_package(catkin_simple REQUIRED)
find_package(catkin REQUIRED COMPONENTS 
    roscpp 
    std_msgs
    nav_msgs
    visualization_msgs 
    geometry_msgs
    message_generation
)

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
find_package(msg_utils REQUIRED)

set(Eigen3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

catkin_package(
    CATKIN_DEPENDS msg_utils
)

include_directories(
    include 
    SYSTEM 
    third_party
    ${catkin_INCLUDE_DIRS} 
    ${Eigen3_INCLUDE_DIRS} 
    ${PCL_INCLUDE_DIRS}
    ${MSG_UTILS_INCLUDE_DIRS}
)

# set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} -O3 -Wall") # -Wextra -Werror


catkin_simple()

cs_add_executable ( random_complex 
    src/random_complex_generator.cpp )

target_link_libraries( random_complex
    ${catkin_LIBRARIES}
    ${PCL_LIBRARIES} )  

cs_install()

cs_export()
