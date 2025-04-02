#ifndef _MPC_SOLVER_H_
#define _MPC_SOLVER_H_

#include "gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <random>
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

#include <boost/math/special_functions.hpp>

typedef struct
{
    Eigen::Vector3d point;
    double a, b, c;
} Obstacle;

class MPCPathFinder
{
private:
    // MPC config
    static constexpr int N = 20;
    double h = 0.05;
    double vsamp = 2.5;
    double g = 9.81;
    double rad = 0.2;
    double delta = 0.03;
    int kapa = 5;

    double posmin[3] = {-15, -15, 0.0};
    double posmax[3] = {15, 15, 5};
    double vmax = 100;
    double amax[3] = {2 * g, 2 * g, g};
    double jmax = 40;
    double Rx = 200, RN = 100, Ru = 0.5, Rj = 0.1, Re = 2.0, Rc = 20.0;
    double thresh_dist = 5.0;
    double eps = 0.2;

    Eigen::MatrixXd A, B;
    Eigen::MatrixXd Q, Gamma, Sigma;
    double qp = 0.05, qv = 0.1;

public:
    MPCPathFinder() {};
    ~MPCPathFinder() {};

    double c = 0.1;
    double w = 0.01;

    int agent_id;
    int agent_num;

    bool flag_ref = false;
    bool flag_init_mpc = false;
    bool flag_failure = false;

    std::vector<Eigen::Vector3d> ref_nodes;

    std::vector<Eigen::MatrixXd> nbr_pred_traj;

    // 记录
    std::vector<Eigen::Vector3d> pos_nodes;
    std::vector<Eigen::Vector3d> v_nodes;
    std::vector<Eigen::Vector3d> posN_nodes;

    std::vector<Obstacle> obs;
    std::vector<Eigen::Vector3d> jps_path;

    double pos[3] = {0, 0, 0};
    double v[3] = {0, 0, 0};
    double a[3] = {0, 0, 0};
    double jerk[3] = {0, 0, 0};
    double goal[3] = {3, 3, 1};
    double pos_n[3];

    double posN[N][3] = {{0}};
    double vN[N][3] = {{0}};
    double aN[N][3] = {{0}};
    double jerkN[N][3] = {{0}};

    int refNum = 0;
    double refposN[N][3];
    int num = 0, iteration = 0, fail = 0;

    ros::Publisher jps_path_pub;
    ros::Publisher es_pub;
    ros::Publisher poly_pub;

    bool check_arrived(double *pos, double *goal);
    bool check_update(double *pos, double *ref);
    bool check_update(double *pos_0, double *ref_0, double *pos_n, double *ref_n);

    // Get Ref
    void GetRefPath(std::vector<Eigen::Vector3d> nodes);
    std::vector<Obstacle> GetNearObs(std::vector<Obstacle> obs, Eigen::Vector3d pos);

    // MPCSolver
    void InitMPCSolver();
    void MPCSolver();
};

#endif