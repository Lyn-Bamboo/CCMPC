#include "gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/PointCloud2.h>
#include <Eigen/Eigen>
#include "Eigen/Dense"

#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>

#include <sensor_msgs/point_cloud_conversion.h>
#include <msg_utils/pred_traj.h>
#include <msg_utils/obsData.h>
#include <msg_utils/obsArray.h>

#include "Astar_searcher.h"
#include "JPS_searcher.h"
#include "backward.hpp"
#include "MPC_solver.h"

class Planner : public JPSPathFinder, public MPCPathFinder
{
private:
    // map param
    double _resolution, _inv_resolution, _cloud_margin;
    double _x_size, _y_size, _z_size;

    Eigen::Vector3d _map_lower, _map_upper;
    int _max_x_id, _max_y_id, _max_z_id;

    double safe_dist = 0.5;

    // ros related
    ros::Subscriber map_sub;
    ros::Subscriber obs_data_sub;
    ros::Subscriber nbr_pred_traj_sub;

    ros::Publisher jps_nodes_pub, visited_nodes_pub;
    ros::Publisher mpc_path_pub, pred_path_pub;
    ros::Publisher pred_traj_share_pub;

public:
    bool flag_map = false;
    bool flag_obs = false;
    // bool flag_ref = false;
    bool flag_arrived = false;

    Eigen::Vector3d start_pt, goal_pt;

    double jps_search_time;
    double ref_path_time;
    std::vector<double> mpc_solving_time;

    // 订阅地图信息
    void MapCallBack(const sensor_msgs::PointCloud2ConstPtr &pointcloud_map);
    // 订阅邻机预测轨迹
    void NBRPredTrajCallBack(const msg_utils::pred_trajConstPtr &msg);
    void ObsDataCallBack(const msg_utils::obsArrayConstPtr &msg);

    // 发布路径话题
    void visJPSNodes(std::vector<Eigen::Vector3d> nodes, Eigen::Vector3d start_pt);
    void visVisitedNode(std::vector<Eigen::Vector3d> nodes);
    void visMPCPath(std::vector<Eigen::Vector3d> nodes);
    void visPredTraj(std::vector<Eigen::Vector3d> nodes);

    // 共享预测轨迹
    void SharePredTraj(std::vector<Eigen::Vector3d> nodes);

    // Find JPS Path
    void JPSPathFinding();

    // Find MPC Path
    void MPCPathFinding();

    Planner(ros::NodeHandle *nh_prevate, ros::NodeHandle *nh);
    Planner(ros::NodeHandle *nh, int id);
    ~Planner() {};
};