cmake_minimum_required(VERSION 2.8.3)
project(planner)


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
find_package(cmake_modules REQUIRED)
find_package(msg_utils REQUIRED)
find_package(Boost REQUIRED COMPONENTS math_c99)



set(Eigen3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

catkin_package(
 INCLUDE_DIRS include
 CATKIN_DEPENDS msg_utils
)

include(./FindGUROBI.cmake)

include_directories(
    include 
    SYSTEM 
    third_party
    ${catkin_INCLUDE_DIRS} 
    ${Eigen3_INCLUDE_DIRS} 
    ${PCL_INCLUDE_DIRS}
    ${GUROBI_INCLUDE_DIRS}
    ${MSG_UTILS_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)


catkin_simple()



cs_add_executable( swarm 
    src/swarm.cpp
    src/Astar_searcher.cpp
    src/JPS_utils.cpp
    src/JPS_searcher.cpp
    src/MPC_solver.cpp
    )

target_link_libraries( swarm 
    ${catkin_LIBRARIES}
    ${PCL_LIBRARIES} 
    ${GUROBI_LIBRARIES}
    ${Boost_LIBRARIES}
)


cs_install()

cs_export()
