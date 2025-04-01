cmake_minimum_required(VERSION 2.8.3)
project(msg_utils)

find_package(catkin REQUIRED COMPONENTS
  #bspline_opt
  #path_searching
  roscpp
  std_msgs
  geometry_msgs
  message_generation
  #cv_bridge
)

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)

# Generate messages in the 'msg' folder
add_message_files(
  FILES
  pred_traj.msg
  obsData.msg
  obsArray.msg
  )
  
# Generate added messages and services with any dependencies listed here
generate_messages(
  DEPENDENCIES
  std_msgs
  geometry_msgs
)

catkin_package(
 CATKIN_DEPENDS message_runtime
)

include_directories( 
    SYSTEM 
    include 
    ${catkin_INCLUDE_DIRS}
    ${Eigen3_INCLUDE_DIRS} 
    ${PCL_INCLUDE_DIRS}
)

link_directories(${PCL_LIBRARY_DIRS})

