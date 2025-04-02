#include <iostream>
#include <fstream>
#include <math.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/PointCloud2.h>

#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>

#include <decomp_ros_utils/data_ros_utils.h>
#include <decomp_util/ellipsoid_decomp.h>

#include "Astar_searcher.h"
#include "JPS_searcher.h"
#include "backward.hpp"
#include "MPC_solver.h"

namespace backward
{
    backward::SignalHandling sh;
}

// map param
double _resolution, _inv_resolution, _cloud_margin;
double _x_size, _y_size, _z_size;

Eigen::Vector3d _map_lower, _map_upper;
int _max_x_id, _max_y_id, _max_z_id;

bool flag_map = false;
bool flag_ref = false;

// start and goal
Eigen::Vector3d start_pt, goal_pt;
// mpc ref path
std::vector<Eigen::Vector3d> ref_nodes;
// SC
vec_Vec3f obs;
vec_Vec3f jps_path;

std::vector<MatDNf<3>> AVec;
std::vector<VecDf> bVec;

// ros related
ros::Subscriber map_sub;
ros::Publisher jps_nodes_vis_pub, visited_nodes_vis_pub;
ros::Publisher jps_path_pub, mpc_path_pub;
ros::Publisher es_pub;
ros::Publisher poly_pub;

JPSPathFinder *_jps_path_finder = new JPSPathFinder();
MPCPathFinder *_mpc_path_finder = new MPCPathFinder();

// Set Obs
void MapCallBack(const sensor_msgs::PointCloud2 &pointcloud_map);

// Find JPS Path
void JPSPathFinding(const Eigen::Vector3d start_pt, const Eigen::Vector3d goal_pt);
void visGridPath(std::vector<Eigen::Vector3d> nodes, Eigen::Vector3d start_pt); // Pub JPS Nodes
void visVisitedNode(std::vector<Eigen::Vector3d> nodes);

// Get Ref Path
void GetRefPath(vec_Vec3f jps_nodes);
// Get SC
void GetSC(visualization_msgs::Marker &markerMsg);

// Find MPC Path
void MPCsolver();
void visMPCPath(std::vector<Eigen::Vector3d> nodes);

void MapCallBack(const sensor_msgs::PointCloud2 &pointcloud_map)
{
    if (flag_map)
        return;

    pcl::PointCloud<pcl::PointXYZ> cloud;
    pcl::PointCloud<pcl::PointXYZ> cloud_vis;
    sensor_msgs::PointCloud2 map_vis;

    pcl::fromROSMsg(pointcloud_map, cloud);

    if ((int)cloud.points.size() == 0)
        return;

    pcl::PointXYZ pt;
    obs.resize(cloud.points.size());
    for (int idx = 0; idx < (int)cloud.points.size(); idx++)
    {
        pt = cloud.points[idx];

        // set obstalces into grid map for path planning
        // 将障碍物信息设置进入栅格化地图，为后续路径规划做准备
        _jps_path_finder->setObs(pt.x, pt.y, pt.z);

        // 为后续求解安全走廊作准备
        obs[idx](0) = pt.x;
        obs[idx](1) = pt.y;
        obs[idx](2) = pt.z;
    }

    flag_map = true;

    JPSPathFinding(start_pt, goal_pt);
}

void JPSPathFinding(const Eigen::Vector3d start_pt, const Eigen::Vector3d goal_pt)
{
    // Call JPS to search for a path
    _jps_path_finder->JPSGraphSearch(start_pt, goal_pt);

    // Retrieve the path
    auto grid_path = _jps_path_finder->getPath();
    auto visited_nodes = _jps_path_finder->getVisitedNodes();

    // Visualize the result
    visGridPath(grid_path, start_pt);
    visVisitedNode(visited_nodes);

    // Reset map for next call
    // _jps_path_finder->resetUsedGrids();
}

void visGridPath(std::vector<Eigen::Vector3d> nodes, Eigen::Vector3d start_pt)
{
    visualization_msgs::Marker node_vis;
    node_vis.header.frame_id = "world";
    node_vis.header.stamp = ros::Time::now();
    node_vis.ns = "demo_node/jps_nodes";

    node_vis.type = visualization_msgs::Marker::CUBE_LIST;
    node_vis.action = visualization_msgs::Marker::ADD;
    node_vis.id = 0;

    node_vis.pose.orientation.x = 0.0;
    node_vis.pose.orientation.y = 0.0;
    node_vis.pose.orientation.z = 0.0;
    node_vis.pose.orientation.w = 1.0;

    node_vis.color.a = 0.5;
    node_vis.color.r = 1.0;
    node_vis.color.g = 0.0;
    node_vis.color.b = 0.0;

    node_vis.scale.x = _resolution;
    node_vis.scale.y = _resolution;
    node_vis.scale.z = _resolution;

    geometry_msgs::Point pt;
    pt.x = start_pt(0);
    pt.y = start_pt(1);
    pt.z = start_pt(2);
    node_vis.points.push_back(pt);

    for (int i = 0; i < int(nodes.size()); i++)
    {
        Eigen::Vector3d coord = nodes[i];
        pt.x = coord(0);
        pt.y = coord(1);
        pt.z = coord(2);

        node_vis.points.push_back(pt);
    }

    jps_nodes_vis_pub.publish(node_vis);
    GetSC(node_vis);
}

void visVisitedNode(std::vector<Eigen::Vector3d> nodes)
{
    visualization_msgs::Marker node_vis;
    node_vis.header.frame_id = "world";
    node_vis.header.stamp = ros::Time::now();
    node_vis.ns = "demo_node/expanded_nodes";
    node_vis.type = visualization_msgs::Marker::CUBE_LIST;
    node_vis.action = visualization_msgs::Marker::ADD;
    node_vis.id = 0;

    node_vis.pose.orientation.x = 0.0;
    node_vis.pose.orientation.y = 0.0;
    node_vis.pose.orientation.z = 0.0;
    node_vis.pose.orientation.w = 1.0;
    node_vis.color.a = 0.5;
    node_vis.color.r = 0.0;
    node_vis.color.g = 0.0;
    node_vis.color.b = 1.0;

    node_vis.scale.x = _resolution;
    node_vis.scale.y = _resolution;
    node_vis.scale.z = _resolution;

    geometry_msgs::Point pt;
    for (int i = 0; i < int(nodes.size()); i++)
    {
        Eigen::Vector3d coord = nodes[i];
        pt.x = coord(0);
        pt.y = coord(1);
        pt.z = coord(2);

        node_vis.points.push_back(pt);
    }

    visited_nodes_vis_pub.publish(node_vis);
}

void GetRefPath(vec_Vec3f jps_nodes)
{
    ref_nodes.clear();
    Eigen::Vector3d start, end, coord;

    for (int i = 0; i < int(jps_nodes.size() - 1); ++i)
    {
        start = jps_nodes[i];
        end = jps_nodes[i + 1];
        double dx = end(0) - start(0);
        double dy = end(1) - start(1);
        double dz = end(2) - start(2);
        double segment_length = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2));

        // 计算分段上等间距采样点数量
        int num_samples = std::max(1, static_cast<int>(segment_length / (vsamp * h)));

        // 计算分段上等间距采样点
        for (int j = 0; j <= num_samples; ++j)
        {
            coord(0) = start(0) + dx * j / num_samples;
            coord(1) = start(1) + dy * j / num_samples;
            coord(2) = start(2) + dz * j / num_samples;
            ref_nodes.push_back(coord);
        }
    }
    ref_nodes.push_back(jps_nodes.back());
}

void GetSC(visualization_msgs::Marker &markerMsg)
{
    // Get JPS
    jps_path.resize(markerMsg->points.size());
    for (int i = 0; i < int(markerMsg->points.size()); i++)
    {
        jps_path[i](0) = markerMsg->points[i].x;
        jps_path[i](1) = markerMsg->points[i].y;
        jps_path[i](2) = markerMsg->points[i].z;
    }
    ROS_INFO("NUMBER OF JPS POINTS IS %ld", jps_path.size());

    nav_msgs::Path jps_path_msg = DecompROS::vec_to_path(jps_path);
    jps_path_msg.header.frame_id = "world";
    jps_path_pub.publish(jps_path_msg);

    // Get Ref
    GetRefPath(jps_path);

    // Get SC
    EllipsoidDecomp3D decomp_util;
    decomp_util.set_obs(obs);
    decomp_util.set_local_bbox(Vec3f(1, 2, 1));
    decomp_util.dilate(jps_path); // Set max iteration number of 10, do fix the path

    // Publish visualization msgs
    decomp_ros_msgs::EllipsoidArray es_msg = DecompROS::ellipsoid_array_to_ros(decomp_util.get_ellipsoids());
    es_msg.header.frame_id = "world";
    es_pub.publish(es_msg);

    decomp_ros_msgs::PolyhedronArray poly_msg = DecompROS::polyhedron_array_to_ros(decomp_util.get_polyhedrons());
    poly_msg.header.frame_id = "world";
    poly_pub.publish(poly_msg);

    // Convert to inequality constraints Ax < b
    AVec.clear();
    bVec.clear();
    auto polys = decomp_util.get_polyhedrons();
    for (size_t i = 0; i < jps_path.size() - 1; i++)
    {
        const auto pt_inside = (jps_path[i] + jps_path[i + 1]) / 2;
        LinearConstraint3D cs(pt_inside, polys[i].hyperplanes());
        AVec.push_back(cs.A());
        bVec.push_back(cs.b());
    }
    ROS_INFO("AVec SIZE IS %ld", AVec.size());

    flag_ref = true;
}

void MPCsolver()
{
}

void visMPCPath(std::vector<Eigen::Vector3d> nodes)
{

    visualization_msgs::Marker node_vis;
    node_vis.header.frame_id = "world";
    node_vis.header.stamp = ros::Time::now();
    node_vis.ns = "mpc/mpc_path_nodes";
    // node_vis.type = visualization_msgs::Marker::SPHERE;
    node_vis.type = visualization_msgs::Marker::LINE_STRIP;
    node_vis.action = visualization_msgs::Marker::ADD;
    node_vis.id = 0;

    node_vis.pose.orientation.x = 0.0;
    node_vis.pose.orientation.y = 0.0;
    node_vis.pose.orientation.z = 0.0;
    node_vis.pose.orientation.w = 1.0;

    node_vis.color.a = 1.0;
    node_vis.color.r = 0.0;
    node_vis.color.g = 0.0;
    node_vis.color.b = 0.0;

    node_vis.scale.x = 0.1;
    // node_vis.scale.y = 0.1;
    // node_vis.scale.z = 0.1;

    geometry_msgs::Point pt;
    for (int i = 0; i < int(nodes.size()); i++)
    {
        Eigen::Vector3d coord = nodes[i];
        pt.x = coord(0);
        pt.y = coord(1);
        pt.z = coord(2);

        // std::cout << "i = " << i << ": x:" << pt.x << "  y:" << pt.y << " z:" << pt.z << std::endl;

        node_vis.points.push_back(pt);
    }

    mpc_path_pub.publish(node_vis);
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "demo_node");
    ros::NodeHandle nh("~");

    // 订阅地图信息的回调函数
    map_sub = nh.subscribe("/random_complex/global_map", 1, MapCallBack);

    jps_nodes_vis_pub = nh.advertise<visualization_msgs::Marker>("jps_nodes_vis", 1);
    visited_nodes_vis_pub = nh.advertise<visualization_msgs::Marker>("visited_nodes_vis", 1);

    es_pub = nh.advertise<decomp_ros_msgs::EllipsoidArray>("ellipsoid_array", 1, true);
    poly_pub = nh.advertise<decomp_ros_msgs::PolyhedronArray>("polyhedron_array", 1, true);

    jps_path_pub = nh.advertise<nav_msgs::Path>("jps_path_vis", 1, true);
    mpc_path_pub = nh.advertise<visualization_msgs::Marker>("mpc_path_vis", 1);

    nh.param("map/cloud_margin", _cloud_margin, 0.0);
    nh.param("map/resolution", _resolution, 0.2);

    nh.param("map/x_size", _x_size, 50.0);
    nh.param("map/y_size", _y_size, 50.0);
    nh.param("map/z_size", _z_size, 5.0);

    nh.param("planning/start_x", start_pt(0), 0.0);
    nh.param("planning/start_y", start_pt(1), 0.0);
    nh.param("planning/start_z", start_pt(2), 0.0);

    nh.param("planning/goal_x", goal_pt(0), 0.0);
    nh.param("planning/goal_y", goal_pt(1), 0.0);
    nh.param("planning/goal_z", goal_pt(2), 0.0);

    _map_lower << -_x_size / 2.0, -_y_size / 2.0, 0.0;
    _map_upper << +_x_size / 2.0, +_y_size / 2.0, _z_size;

    _inv_resolution = 1.0 / _resolution;

    _max_x_id = (int)(_x_size * _inv_resolution);
    _max_y_id = (int)(_y_size * _inv_resolution);
    _max_z_id = (int)(_z_size * _inv_resolution);

    // 定义了结构体 JPSPathFinder 变量_jps_path_finder
    // 该结构体存储、实现了 JPS 路径规划所需的所有信息和功能
    _jps_path_finder = new JPSPathFinder();
    _jps_path_finder->initGridMap(_resolution, _map_lower, _map_upper, _max_x_id, _max_y_id, _max_z_id);

    ros::Rate rate(100);
    bool status = ros::ok();
    while (status && flag_ref)
    {

        ros::spinOnce();
        status = ros::ok();
        rate.sleep();
    }

    delete _jps_path_finder;
    return 0;
}
