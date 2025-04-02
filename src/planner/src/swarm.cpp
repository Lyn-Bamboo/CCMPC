#include "planner.h"

namespace backward
{
    backward::SignalHandling sh;
}
double _x_size;
void Planner::MapCallBack(const sensor_msgs::PointCloud2ConstPtr &pointcloud_map)
{
    if (flag_map)
        return;

    pcl::PointCloud<pcl::PointXYZ> cloud;
    pcl::PointCloud<pcl::PointXYZ> cloud_vis;
    sensor_msgs::PointCloud2 map_vis;

    pcl::fromROSMsg(*pointcloud_map, cloud);

    if ((int)cloud.points.size() == 0)
        return;

    pcl::PointXYZ pt;
    for (int idx = 0; idx < (int)cloud.points.size(); idx++)
    {
        pt = cloud.points[idx];

        // set obstalces into grid map for path planning
        // 将障碍物信息设置进入栅格化地图，为后续JPS路径规划做准备
        setObs(pt.x, pt.y, pt.z);
    }

    flag_map = true;
    ROS_INFO("MAP--------------------------------------------------------");
}

void Planner::ObsDataCallBack(const msg_utils::obsArrayConstPtr &msg)
{
    if (flag_obs)
        return;

    obs.resize(msg->obsArray.size());
    for (int i = 0; i < (int)msg->obsArray.size(); i++)
    {
        obs[i].point(0) = msg->obsArray[i].x;
        obs[i].point(1) = msg->obsArray[i].y;
        obs[i].point(2) = msg->obsArray[i].z;
        obs[i].a = msg->obsArray[i].a;
        obs[i].b = msg->obsArray[i].b;
        obs[i].c = msg->obsArray[i].c;
    }
    ROS_INFO("THERE ARE %ld OBS", msg->obsArray.size());

    flag_obs = true;
}

void Planner::NBRPredTrajCallBack(const msg_utils::pred_trajConstPtr &msg)
{
    int id = msg->agent_id;

    // 检查接收到的消息是否是来自本机的，如果是，则直接返回，不进行处理
    if (id == agent_id)
        return;

    // 检查接收到的消息代表的起点与本机当前位置之间的距离是否超过了规定的阈值
    bool safe = true;
    for (int i = 0; i < msg->points.size(); i++)
    {
        Eigen::Vector3d pj(msg->points[i].x, msg->points[i].y, msg->points[i].z);
        if ((pj - posN_nodes[i]).norm() <= 6 * safe_dist)
        {
            // 存储
            nbr_pred_traj[id].resize(msg->points.size(), 3);
            for (int j = 0; j < msg->points.size(); j++)
            {
                nbr_pred_traj[id](j, 0) = msg->points[j].x;
                nbr_pred_traj[id](j, 1) = msg->points[j].y;
                nbr_pred_traj[id](j, 2) = msg->points[j].z;
            }
            ROS_INFO("AGENT %d GET NBR PRED TRAJECTORY!!!!!!!!!!!!!!!!", agent_id);
            safe = false;
            return;
        }
    }
    if (safe)
        nbr_pred_traj[id].resize(0, 0); // if the current agent is too far to the received agent.
}

void Planner::visJPSNodes(std::vector<Eigen::Vector3d> nodes, Eigen::Vector3d start_pt)
{
    visualization_msgs::Marker node_vis;
    node_vis.header.frame_id = "world";
    node_vis.header.stamp = ros::Time::now();
    node_vis.ns = "/jps_nodes";

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

    jps_nodes_pub.publish(node_vis);
}

void Planner::visVisitedNode(std::vector<Eigen::Vector3d> nodes)
{
    visualization_msgs::Marker node_vis;
    node_vis.header.frame_id = "world";
    node_vis.header.stamp = ros::Time::now();
    node_vis.ns = "/expanded_nodes";
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

    visited_nodes_pub.publish(node_vis);
}

void Planner::visMPCPath(std::vector<Eigen::Vector3d> nodes)
{
    nav_msgs::Path node_vis;

    node_vis.header.frame_id = "world";
    node_vis.header.stamp = ros::Time::now();

    geometry_msgs::PoseStamped pt;
    for (int i = 0; i < int(nodes.size()); i++)
    {
        Eigen::Vector3d coord = nodes[i];

        pt.pose.position.x = coord(0);
        pt.pose.position.y = coord(1);
        pt.pose.position.z = coord(2);

        pt.pose.orientation.x = 0.0;
        pt.pose.orientation.y = 0.0;
        pt.pose.orientation.z = 0.0;
        pt.pose.orientation.w = 1.0;

        // std::cout << "i = " << i << ": x:" << pt.x << "  y:" << pt.y << " z:" << pt.z << std::endl;

        node_vis.poses.push_back(pt);
    }

    mpc_path_pub.publish(node_vis);
}

void Planner::visPredTraj(std::vector<Eigen::Vector3d> nodes)
{
    visualization_msgs::Marker node_vis;

    node_vis.header.frame_id = "world";
    node_vis.header.stamp = ros::Time::now();
    node_vis.ns = "/pred_traj";
    node_vis.type = visualization_msgs::Marker::LINE_STRIP;
    node_vis.action = visualization_msgs::Marker::ADD;
    node_vis.id = 0;

    node_vis.pose.orientation.x = 0.0;
    node_vis.pose.orientation.y = 0.0;
    node_vis.pose.orientation.z = 0.0;
    node_vis.pose.orientation.w = 1.0;

    node_vis.color.a = 1.0;
    node_vis.color.r = 0.8;
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

    pred_path_pub.publish(node_vis);
}

void Planner::SharePredTraj(std::vector<Eigen::Vector3d> nodes)
{
    msg_utils::pred_traj traj;

    traj.agent_id = agent_id;

    geometry_msgs::Point pt;
    for (int i = 0; i < int(nodes.size()); i++)
    {
        Eigen::Vector3d coord = nodes[i];
        pt.x = coord(0);
        pt.y = coord(1);
        pt.z = coord(2);

        // std::cout << "i = " << i << ": x:" << pt.x << "  y:" << pt.y << " z:" << pt.z << std::endl;

        traj.points.push_back(pt);
    }

    traj.start_time = ros::Time::now();

    pred_traj_share_pub.publish(traj);
}

void Planner::JPSPathFinding()
{
    // Call JPS to search for a path
    ros::Time t1 = ros::Time::now();

    JPSGraphSearch(start_pt, goal_pt);

    ros::Time t2 = ros::Time::now();
    jps_search_time = (t2 - t1).toSec() * 1000.0;

    // Retrieve the path
    auto grid_path = getPath();
    auto visited_nodes = getVisitedNodes();

    // Visualize the result
    visJPSNodes(grid_path, start_pt);
    visVisitedNode(visited_nodes);

    t1 = ros::Time::now();

    grid_path.insert(grid_path.begin(), start_pt);
    GetRefPath(grid_path);

    t2 = ros::Time::now();
    ref_path_time = (t2 - t1).toSec() * 1000.0;
}

void Planner::MPCPathFinding()
{
    if (!flag_init_mpc)
    {
        InitMPCSolver();
        ROS_INFO("AGENT %d MPC SOLVER INIT !!!", agent_id);
    }

    ROS_INFO("AGENT %d MPC START !!!", agent_id);

    ros::Time t1 = ros::Time::now();

    MPCSolver();

    ros::Time t2 = ros::Time::now();
    mpc_solving_time.push_back((t2 - t1).toSec() * 1000.0);

    ROS_INFO("AGENT %d MPC END !!!", agent_id);

    /************************************* 调试 ********************************************/
    // std::cout << "ref_pos:" << refNum << std::endl;
    // for (int i = 0; i < 3; i++)
    // {
    //     for (int j = 0; j < refNum; j++)
    //     {
    //         std::cout << ref_nodes[j](i) << "  ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // std::cout << "mpc_pos:" << int(pos_nodes.size()) << std::endl;
    // for (int i = 0; i < 3; i++)
    // {
    //     for (int j = 0; j < int(pos_nodes.size()); j++)
    //     {
    //         std::cout << pos_nodes[j](i) << "  ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // std::cout << "mpc_v:" << int(v_nodes.size()) << std::endl;
    // for (int i = 0; i < 3; i++)
    // {
    //     for (int j = 0; j < int(v_nodes.size()); j++)
    //     {
    //         std::cout << v_nodes[j](i) << "  ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // std::cout << "pred_pos:" << int(posN_nodes.size()) << std::endl;
    // for (int i = 0; i < 3; i++)
    // {
    //     for (int j = 0; j < int(posN_nodes.size()); j++)
    //     {
    //         std::cout << posN_nodes[j](i) << "  ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    /************************************* 调试 ********************************************/

    if (check_arrived(pos, goal))
    {
        ROS_INFO("Agent %d Arrive Goal Successfully !!!", agent_id);
        flag_arrived = true;
    }
    if (flag_failure)
    {
        ROS_INFO("Agent %d Planning Failure !!!", agent_id);
    }

    visMPCPath(pos_nodes);
    visPredTraj(posN_nodes);
    SharePredTraj(posN_nodes);
}

Planner::Planner(ros::NodeHandle *nh, int id)
{
    agent_id = id;

    // 订阅地图信息的回调函数
    map_sub = nh->subscribe<sensor_msgs::PointCloud2>("/random_complex/global_map", 1, std::bind(&Planner::MapCallBack, this, std::placeholders::_1));
    obs_data_sub = nh->subscribe<msg_utils::obsArray>("/random_complex/obs_data", 1, std::bind(&Planner::ObsDataCallBack, this, std::placeholders::_1));
    // 订阅邻机预测轨迹
    nbr_pred_traj_sub = nh->subscribe<msg_utils::pred_traj>("/swarm/pred_traj", 20, std::bind(&Planner::NBRPredTrajCallBack, this, std::placeholders::_1));

    // 发布规划路径用于rviz显示
    std::string str = "/agent_" + std::to_string(agent_id) + "_planner/jps_nodes_vis";
    jps_nodes_pub = nh->advertise<visualization_msgs::Marker>(str.c_str(), 10);

    str = "/agent_" + std::to_string(agent_id) + "_planner/visited_nodes_vis";
    visited_nodes_pub = nh->advertise<visualization_msgs::Marker>(str.c_str(), 10);

    str = "/agent_" + std::to_string(agent_id) + "_planner/jps_path_vis";
    jps_path_pub = nh->advertise<nav_msgs::Path>(str.c_str(), 10, true);

    str = "/agent_" + std::to_string(agent_id) + "_planner/mpc_path_vis";
    mpc_path_pub = nh->advertise<nav_msgs::Path>(str.c_str(), 10);

    str = "/agent_" + std::to_string(agent_id) + "_planner/pred_path_vis";
    pred_path_pub = nh->advertise<visualization_msgs::Marker>(str.c_str(), 10);

    // 全局节点共享预测轨迹
    pred_traj_share_pub = nh->advertise<msg_utils::pred_traj>("pred_traj", 1);

    // get param
    // nh->param("agent_id", agent_id, 0);
    nh->param("agent_num", agent_num, 1);

    nh->param("map/cloud_margin", _cloud_margin, 0.0);
    nh->param("map/resolution", _resolution, 0.2);

    nh->param("map/x_size", _x_size, 50.0);
    nh->param("map/y_size", _y_size, 50.0);
    nh->param("map/z_size", _z_size, 5.0);

    nh->param("c", c, 0.1);
    nh->param("w", w, 0.01);

    nh->param("start_" + std::to_string(agent_id) + "/x", start_pt(0), 6.0);
    nh->param("start_" + std::to_string(agent_id) + "/y", start_pt(1), 0.0);
    nh->param("start_" + std::to_string(agent_id) + "/z", start_pt(2), 1.0);

    nh->param("goal_" + std::to_string(agent_id) + "/x", goal_pt(0), -6.0);
    nh->param("goal_" + std::to_string(agent_id) + "/y", goal_pt(1), 0.0);
    nh->param("goal_" + std::to_string(agent_id) + "/z", goal_pt(2), 1.0);

    ROS_INFO("agent_id: %d; start: %f, %f, %f; goal: %f, %f, %f",
             agent_id, start_pt(0), start_pt(1), start_pt(2), goal_pt(0), goal_pt(1), goal_pt(2));
    ROS_INFO("c: %f; w: %f", c, w);

    _map_lower << -_x_size / 2.0, -_y_size / 2.0, 0.0;
    _map_upper << +_x_size / 2.0, +_y_size / 2.0, _z_size;

    _inv_resolution = 1.0 / _resolution;

    _max_x_id = (int)(_x_size * _inv_resolution);
    _max_y_id = (int)(_y_size * _inv_resolution);
    _max_z_id = (int)(_z_size * _inv_resolution);

    initGridMap(_resolution, _map_lower, _map_upper, _max_x_id, _max_y_id, _max_z_id);

    nbr_pred_traj.resize(agent_num);
    for (int id = 0; id < agent_num; id++)
    {
        nbr_pred_traj[id].resize(0, 0);
    }
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "swarm");
    ros::NodeHandle nh("~");
    // ros::NodeHandle nh;

    int num;
    nh.param("agent_num", num, 1);
    std::vector<Planner *> planners; // 使用vector来存储Planner实例的指针

    for (int i = 0; i < num; ++i)
    {
        Planner *planner = new Planner(&nh, i); // 动态分配Planner实例
        planners.push_back(planner);            // 将指针存储到vector中
    }

    ros::Rate rate(100);
    int count = 0;
    int fail = 0;
    double cmp_t = 0;
    while (ros::ok() && count < num)
    {
        for (int i = 0; i < num; i++)
        {
            if (planners[i]->flag_map && !planners[i]->flag_ref)
            {
                ROS_INFO("JPS START ---------------------------------------");
                planners[i]->JPSPathFinding();
                ROS_INFO("JPS END ---------------------------------------");
                count++;
            }
        }
        ros::spinOnce();
    }
    count = 0;
    while (ros::ok() && (count + fail) < num)
    {
        count = 0;
        fail = 0;
        cmp_t = 0;
        for (int i = 0; i < num; i++)
        {
            if (planners[i]->flag_ref && !planners[i]->flag_arrived && !planners[i]->flag_failure)
            {
                planners[i]->MPCPathFinding();
            }
            else if (planners[i]->flag_arrived)
            {
                double max_val = *std::max_element(planners[i]->mpc_solving_time.begin(), planners[i]->mpc_solving_time.end());
                double min_val = *std::min_element(planners[i]->mpc_solving_time.begin(), planners[i]->mpc_solving_time.end());
                double mean = std::accumulate(planners[i]->mpc_solving_time.begin(), planners[i]->mpc_solving_time.end(), 0.0) / planners[i]->mpc_solving_time.size();

                ROS_WARN("AGENT %d ARRIVED SUCCESSFULLY!", planners[i]->agent_id);
                ROS_WARN("jps search time: %f ms, ref path obtain time: %f ms", planners[i]->jps_search_time, planners[i]->ref_path_time);
                ROS_WARN("mpc solving time: min: %f ms, max: %f ms, mean: %f ms", min_val, max_val, mean);
                ROS_INFO("-------------------------------------------------------------");

                cmp_t += mean / (double)num;
                count++;
            }
            else if (planners[i]->flag_failure)
            {
                ROS_ERROR("AGENT %d PLANNING FAILURE!", planners[i]->agent_id);
                fail++;
            }
        }
        ros::spinOnce();
    }

    /****************************************************************************************************/
    // 数据处理

    std::cout << "MEAN COMP TIME:  " << cmp_t << "ms" << std::endl;

    // trajectory length
    std::vector<double> length;
    for (int i = 0; i < num; i++)
    {
        double l = 0;
        for (int j = 1; j < (int)planners[i]->pos_nodes.size(); j++)
        {
            l += (planners[i]->pos_nodes[j] - planners[i]->pos_nodes[j - 1]).norm();
        }
        length.push_back(l);
    }
    double max_l = *std::max_element(length.begin(), length.end());
    double min_l = *std::min_element(length.begin(), length.end());
    double mean_l = std::accumulate(length.begin(), length.end(), 0.0) / length.size();
    ROS_INFO("-------------------------------------------------------------");
    std::cout << "TRAJECTORY LENGTH:" << std::endl;
    for (int i = 0; i < num; i++)
    {
        std::cout << length[i] << "  ";
    }
    std::cout << std::endl;
    std::cout << "min: " << min_l << "m, max: " << max_l << "m, mean: " << mean_l << "m" << std::endl;

    // flight time
    std::vector<double> time;
    for (int i = 0; i < num; i++)
    {
        double tt = 0.1 * planners[i]->pos_nodes.size();
        time.push_back(tt);
    }
    double max_t = *std::max_element(time.begin(), time.end());
    double min_t = *std::min_element(time.begin(), time.end());
    double mean_t = std::accumulate(time.begin(), time.end(), 0.0) / time.size();
    ROS_INFO("-------------------------------------------------------------");
    std::cout << "FIGHT TIME:" << std::endl;
    for (int i = 0; i < num; i++)
    {
        std::cout << time[i] << "  ";
    }
    std::cout << std::endl;
    std::cout << "min: " << min_t << "s, max: " << max_t << "s, mean: " << mean_t << "s" << std::endl;

    // dij: minimum distance
    int point_num = 1e5;
    for (int i = 0; i < num; i++)
    {
        if (point_num > (int)planners[i]->pos_nodes.size())
        {
            point_num = (int)planners[i]->pos_nodes.size();
        }
    }

    double dij = 1e5;
    double min_dij = 1e5;
    for (int i = 0; i < num; i++)
    {
        for (int k = 0; k < point_num; k++)
        {
            dij = 1e5;
            for (int j = 0; j < num; j++)
            {
                if (j == i)
                    continue;
                double dd = (planners[i]->pos_nodes[k] - planners[j]->pos_nodes[k]).norm();
                // std::cout << "ddddij:" << dij << std::endl;
                if (dd < dij)
                    dij = dd;
            }
            if (dij <= 5.0)
            {
                if (dij < min_dij)
                {
                    min_dij = dij;
                }
            }
        }
    }

    std::cout << "MIN D_ij:  " << min_dij << "m" << std::endl;

    // 清空动态分配的内存
    for (auto planner : planners)
    {
        delete planner; // 释放Planner实例的内存
    }

    // 清空vector
    planners.clear(); // 清空vector中的所有元素

    ROS_INFO("-- FINISHED --");
    ros::spin();

    return 0;
}
