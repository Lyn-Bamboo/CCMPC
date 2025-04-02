#include "MPC_solver.h"

bool MPCPathFinder::check_arrived(double *pos, double *goal)
{
    double dist = std::sqrt(std::pow(pos[0] - goal[0], 2) + std::pow(pos[1] - goal[1], 2) + std::pow(pos[2] - goal[2], 2));
    // std::cout << "agent " << agent_id << " goal_dist = " << dist << std::endl;
    return (dist < eps);
}

bool MPCPathFinder::check_update(double *pos, double *ref)
{
    double dist = std::sqrt(std::pow(pos[0] - ref[0], 2) + std::pow(pos[1] - ref[1], 2) + std::pow(pos[2] - ref[2], 2));
    // ROS_INFO("agent %d update_dist = %f", agent_id, dist);
    return (dist < thresh_dist);
}
bool MPCPathFinder::check_update(double *pos_0, double *ref_0, double *pos_n, double *ref_n)
{
    double dist1 = std::sqrt(std::pow(pos_0[0] - ref_0[0], 2) + std::pow(pos_0[1] - ref_0[1], 2) + std::pow(pos_0[2] - ref_0[2], 2));
    double dist2 = std::sqrt(std::pow(pos_n[0] - ref_n[0], 2) + std::pow(pos_n[1] - ref_n[1], 2) + std::pow(pos_n[2] - ref_n[2], 2));
    // ROS_INFO("agent %d update_dist = %f, %f", agent_id, dist1, dist2);
    return (dist1 < thresh_dist || dist2 < thresh_dist);
}

void MPCPathFinder::GetRefPath(std::vector<Eigen::Vector3d> jps_nodes)
{
    ROS_INFO("AGENT %d NUMBER OF JPS POINTS IS %ld", agent_id, jps_nodes.size());

    nav_msgs::Path jps_path_msg;
    for (int i = 0; i < int(jps_nodes.size()); i++)
    {
        geometry_msgs::PoseStamped pose;
        pose.pose.position.x = jps_nodes[i](0);
        pose.pose.position.y = jps_nodes[i](1);
        pose.pose.position.z = jps_nodes[i](2);
        pose.pose.orientation.w = 1.0;
        pose.pose.orientation.x = 0.0;
        pose.pose.orientation.y = 0.0;
        pose.pose.orientation.z = 0.0;

        jps_path_msg.poses.push_back(pose);
    }
    jps_path_msg.header.frame_id = "world";
    jps_path_pub.publish(jps_path_msg);

    Eigen::Vector3d start, end, coord;
    ref_nodes.clear();

    for (int i = 0; i < int(jps_nodes.size() - 1); ++i)
    {
        start = jps_nodes[i];
        end = jps_nodes[i + 1];
        double dx = end(0) - start(0);
        double dy = end(1) - start(1);
        double dz = end(2) - start(2);
        double segment_length = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2));

        // 计算分段上等间距采样点数量
        // int num_samples = std::max(1, static_cast<int>(segment_length / (vsamp * h)));
        double num_samples = std::max(1.0, segment_length / (vsamp * h));
        // dd = num_samples - static_cast<int>(num_samples);+

        // 计算分段上等间距采样点
        for (int j = 0; j < num_samples; ++j)
        {
            coord(0) = start(0) + dx * j / num_samples;
            coord(1) = start(1) + dy * j / num_samples;
            coord(2) = start(2) + dz * j / num_samples;
            ref_nodes.push_back(coord);
        }
    }

    ref_nodes.push_back(jps_nodes.back());

    std::cout << "agent " << agent_id << " jps_pos: " << jps_nodes.size() << std::endl;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < jps_nodes.size(); j++)
        {
            std::cout << jps_nodes[j](i) << "  ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "agent " << agent_id << " ref_pos:" << ref_nodes.size() << std::endl;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < ref_nodes.size(); j++)
        {
            std::cout << ref_nodes[j](i) << "  ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    flag_ref = true;
}

std::vector<Obstacle> MPCPathFinder::GetNearObs(std::vector<Obstacle> obs, Eigen::Vector3d pos)
{
    std::vector<Obstacle> near_obs;
    for (int i = 0; i < (int)obs.size(); i++)
    {
        if ((obs[i].point - pos).norm() <= 3.0)
        {
            near_obs.push_back(obs[i]);
        }
    }
    // ROS_INFO("!!!!!!OBS NUMBER IS %ld", obs.size());
    // ROS_INFO("!!!!!!OBS NUMBER IS %ld", near_obs.size());
    return near_obs;
}

void MPCPathFinder::InitMPCSolver()
{
    // mpc参数
    A.resize(6, 6);
    A << 1, 0, 0, h, 0, 0,
        0, 1, 0, 0, h, 0,
        0, 0, 1, 0, 0, h,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1;
    B.resize(3, 6);
    B << h * h / 2, 0, 0,
        0, h * h / 2, 0,
        0, 0, h * h / 2,
        h, 0, 0,
        0, h, 0,
        0, 0, h;
    Gamma = Q = Eigen::MatrixXd::Zero(6, 6);
    Gamma.diagonal() << qp * qp, qp * qp, qp * qp, qv * qv, qv * qv, qv * qv;
    Q.diagonal() << qp * qp, qp * qp, qp * qp, qv * qv, qv * qv, qv * qv;
    Sigma.resize(3, 3);
    Sigma = Gamma.block(0, 0, 3, 3);

    // 起止点、参考点信息
    pos_nodes.clear();
    v_nodes.clear();
    posN_nodes.resize(N);

    pos_nodes.push_back(ref_nodes[0]);

    refNum = ref_nodes.size();

    for (int i = 0; i < N; i++)
    {
        refposN[i][0] = ref_nodes[i](0);
        refposN[i][1] = ref_nodes[i](1);
        refposN[i][2] = ref_nodes[i](2);

        posN[i][0] = ref_nodes[i](0);
        posN[i][1] = ref_nodes[i](1);
        posN[i][2] = ref_nodes[i](2);
    }

    pos[0] = ref_nodes[0](0);
    pos[1] = ref_nodes[0](1);
    pos[2] = ref_nodes[0](2);

    goal[0] = ref_nodes[refNum - 1](0);
    goal[1] = ref_nodes[refNum - 1](1);
    goal[2] = ref_nodes[refNum - 1](2);

    flag_init_mpc = true;
}

void MPCPathFinder::MPCSolver()
{
    if (!check_arrived(pos, goal) && fail <= 1 * N)
    {

        // std::cout << "---------------------------------------------------------------------------------------" << std::endl;
        // ROS_INFO("NUM = %d, ITERATION = %d", num, iteration);

        // Solve MPC Problem
        try
        {
            int cc = 0;

            // Create an environment
            GRBEnv env = GRBEnv(true);

            env.set("LogFile", "mpc.log");
            env.start();

            // Create an empty model
            GRBModel model = GRBModel(env);
            // 设置线程数
            model.set(GRB_IntParam_Threads, 1);
            // 启用输出信息
            model.getEnv().set(GRB_IntParam_OutputFlag, 1);
            // 设置求解时间
            // model.getEnv().set(GRB_DoubleParam_TimeLimit, 0.1);

            // 设置求解重点：0 默认，均衡搜寻可行解和证明最优；1 侧重快速找到可行解； 2 侧重证明最优； 3 侧重界的提升（发现界提升缓慢）
            // GRB_MIPFOCUS_BALANCED：默认设置，求解器在寻找最优解和快速求解之间寻求平衡。
            // GRB_MIPFOCUS_FEASIBILITY：求解器侧重于快速找到一个可行解，而不是最优解。
            // GRB_MIPFOCUS_OPTIMALITY：求解器侧重于找到最优解，可能会花费更长的时间。
            // GRB_MIPFOCUS_BESTBOUND：求解器侧重于改进目标函数的界限（即，找到一个更好的边界值）。
            // GRB_MIPFOCUS_HIDDENFEAS：求解器侧重于寻找隐藏的可行解，这通常在特定情况下使用。
            model.getEnv().set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_FEASIBILITY);

            // Create variables
            GRBVar posNvar[N][3];
            GRBVar vNvar[N][3];
            GRBVar aNvar[N][3];
            GRBVar jerkNvar[N][3];

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    posNvar[i][j] = model.addVar(posmin[j], posmax[j], 0.0, GRB_CONTINUOUS, "posNvar_" + std::to_string(i) + "_" + std::to_string(j));
                    vNvar[i][j] = model.addVar(-vmax, vmax, 0.0, GRB_CONTINUOUS, "vNvar_" + std::to_string(i) + "_" + std::to_string(j));
                    aNvar[i][j] = model.addVar(-amax[j], amax[j], 0.0, GRB_CONTINUOUS, "aNvar_" + std::to_string(i) + "_" + std::to_string(j));
                    jerkNvar[i][j] = model.addVar(-jmax, jmax, 0.0, GRB_CONTINUOUS, "jerkNvar_" + std::to_string(i) + "_" + std::to_string(j));
                }
            }

            // set start value
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    posNvar[i][j].set(GRB_DoubleAttr_Start, refposN[i][j]);
                }
            }

            // Add Constrains 0

            for (int j = 0; j < 3; j++)
            {
                model.addConstr(vNvar[N - 1][j] == 0, "c" + std::to_string(cc++));
                model.addConstr(aNvar[N - 1][j] == 0, "c" + std::to_string(cc++));
            }

            // Add Constrains 1
            // Xk+1 = A*Xk + B
            // N = 1

            for (int j = 0; j < 3; j++)
            {
                model.addConstr(posNvar[0][j] == pos[j] + h * v[j] + h * h / 2 * a[j], "c" + std::to_string(cc++));
                model.addConstr(vNvar[0][j] == v[j] + h * a[j], "c" + std::to_string(cc++));
            }
            // N > 1
            for (int i = 1; i < N; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    model.addConstr(posNvar[i][j] == posNvar[i - 1][j] + h * vNvar[i - 1][j] + h * h / 2 * aNvar[i - 1][j], "c" + std::to_string(cc++));
                    model.addConstr(vNvar[i][j] == vNvar[i - 1][j] + h * aNvar[i - 1][j], "c" + std::to_string(cc++));
                }
            }

            // Add Constrains 2
            // ||pi-po|| < b
            Eigen::MatrixXd Gamma1(6, 6), Gamma2(6, 6);
            Eigen::MatrixXd Sigma1(3, 3), Sigma2(3, 3);
            Gamma1 = Gamma;

            std::vector<Obstacle> near_obs = GetNearObs(obs, pos_nodes.back());
            int obs_num = (int)near_obs.size();
            ROS_INFO("AGENT %d : OBS NUMBER IS %d", agent_id, obs_num);

            for (int i = 0; i < N; i++)
            {
                if (!num)
                    break;
                Gamma2 = A * Gamma1 * A.transpose() + Q;
                Sigma2 = Gamma2.block(0, 0, 3, 3);

                if (obs_num)
                {
                    for (int idx = 0; idx < obs_num; idx++)
                    {
                        // Omega
                        Eigen::MatrixXd Omega = Eigen::MatrixXd::Zero(3, 3);
                        Omega.diagonal() << 1.0 / (near_obs[idx].a + rad), 1.0 / (near_obs[idx].b + rad), 1.0 / (near_obs[idx].c + rad);
                        // a_io
                        Eigen::Vector3d a_io = Omega * (posN_nodes[i] - near_obs[idx].point) / (Omega * (posN_nodes[i] - near_obs[idx].point)).norm();
                        // b_io
                        Eigen::MatrixXd tt = 2 * a_io.transpose() * Omega * Sigma2 * Omega.transpose() * a_io;
                        double b_io = 1 + boost::math::erf_inv(1 - 2 * delta) * std::sqrt(tt(0, 0));
                        // std::cout << "b_io: " << b_io << std::endl;
                        // a_io*(pi-po) >= b_io
                        model.addConstr(a_io(0) * Omega(0, 0) * (posNvar[i][0] - near_obs[idx].point(0)) +
                                                a_io(1) * Omega(1, 1) * (posNvar[i][1] - near_obs[idx].point(1)) +
                                                a_io(2) * Omega(2, 2) * (posNvar[i][2] - near_obs[idx].point(2)) >=
                                            b_io,
                                        "c" + std::to_string(cc++));
                    }
                }
                Gamma1 = Gamma2;
            }

            // Add Constrains 3
            // ||pi-pj|| < ri + rj (b_ij)

            int cnt = 0;
            for (int k = 0; k < agent_num; k++)
            {
                if (!num)
                    break;

                if (k == agent_id)
                    continue;

                if (nbr_pred_traj[k].rows() == 0 || nbr_pred_traj[k].cols() == 0)
                    continue;

                Gamma1 = Gamma;
                for (int i = 0; i < N; i++)
                {
                    Gamma2 = A * Gamma1 * A.transpose() + Q;
                    Sigma2 = Gamma2.block(0, 0, 3, 3);
                    // a_ij
                    Eigen::Vector3d pj = nbr_pred_traj[k].row(i);
                    Eigen::Vector3d a_ij = (posN_nodes[i] - pj) / (posN_nodes[i] - pj).norm();
                    // b_ij
                    Eigen::MatrixXd tt = 2 * a_ij.transpose() * (Sigma2 + Sigma2) * a_ij;
                    double b_ij = 2 * rad + boost::math::erf_inv(1 - 2 * delta) * std::sqrt(tt(0, 0));
                    // std::cout << "b_ij: " << b_ij << std::endl;
                    // a_io*(pi-po) >= b_io
                    model.addConstr(a_ij(0) * (posNvar[i][0] - pj(0)) +
                                            a_ij(1) * (posNvar[i][1] - pj(1)) +
                                            a_ij(2) * (posNvar[i][2] - pj(2)) >=
                                        b_ij,
                                    "c" + std::to_string(cc++));

                    Gamma1 = Gamma2;
                }
                cnt++;
            }

            ROS_INFO("AGENT %d : NBR AGENT NUMBER IS %d", agent_id, cnt);

            //
            //
            // Set objective: obj = Jx + Ju + JN;
            // 设置目标函数
            GRBQuadExpr obj = 0;
            GRBQuadExpr Jx = 0;
            GRBQuadExpr JN = 0;
            GRBQuadExpr Ju = 0;
            GRBQuadExpr Jj = 0;

            GRBQuadExpr Je = 0;
            GRBQuadExpr Jc = 0;

            // JN = ||pos_n-ref_n||RN
            for (int j = 0; j < 3; j++)
            {
                JN += (posNvar[N - 1][j] - refposN[N - 1][j]) * RN * (posNvar[N - 1][j] - refposN[N - 1][j]);
            }
            // Jx = ||pos-ref||Rx
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    Jx += (posNvar[i][j] - refposN[i][j]) * Rx * (posNvar[i][j] - refposN[i][j]);
                }
            }

            // // Je = ||posN(N-kapa:N)-goal||Re
            // double ee = (pos[0] - goal[0]) * (pos[0] - goal[0]) + (pos[1] - goal[1]) * (pos[1] - goal[1]) + (pos[2] - goal[2]) * (pos[2] - goal[2]);
            // for (int i = 1; i <= kapa; i++)
            // {
            //     for (int j = 0; j < 3; j++)
            //     {
            //         Je += (posNvar[N - i][j] - goal[j]) * Re * (posNvar[N - i][j] - goal[j]);
            //     }
            // }

            // Ju = ||a||Ru
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    Ju += aNvar[i][j] * Ru * aNvar[i][j];
                }
            }
            // Jj = ||a(k+1)-a(k)||Rj
            for (int i = 0; i < N - 1; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    if (!i)
                    {
                        Jj += (aNvar[i][j] - a[j]) * Rj * (aNvar[i][j] - a[j]);
                    }
                    Jj += (aNvar[i + 1][j] - aNvar[i][j]) * Rj * (aNvar[i + 1][j] - aNvar[i][j]);
                }
            }

            // 设置为简化的目标函数
            obj = Jx + Ju + JN + Jj + Je + Jc;
            model.setObjective(obj, GRB_MINIMIZE);

            model.write("model.lp");

            //
            //
            // Optimize model
            model.optimize();

            // 检查模型状态
            int status = model.get(GRB_IntAttr_Status);

            switch (status)
            {
            case GRB_OPTIMAL:
                std::cout << "Optimal solution found." << std::endl;
                // ... 可以访问解决方案属性，如 ObjVal, XVar 等 ...

                /*********************************** no reference *****************************************/
                // fail = 0;
                // iteration++;

                // for (int j = 0; j < 3; j++)
                // {
                //     pos[j] = posNvar[0][j].get(GRB_DoubleAttr_X);
                //     v[j] = vNvar[0][j].get(GRB_DoubleAttr_X);
                //     a[j] = aNvar[0][j].get(GRB_DoubleAttr_X);
                //     jerk[j] = jerkNvar[0][j].get(GRB_DoubleAttr_X);
                //     pos_n[j] = posNvar[N - 1][j].get(GRB_DoubleAttr_X);
                // }

                // // update xN
                // for (int i = 0; i < N; i++)
                // {
                //     for (int j = 0; j < 3; j++)
                //     {
                //         posN[i][j] = posNvar[i][j].get(GRB_DoubleAttr_X);
                //         vN[i][j] = vNvar[i][j].get(GRB_DoubleAttr_X);
                //         aN[i][j] = aNvar[i][j].get(GRB_DoubleAttr_X);
                //         jerkN[i][j] = jerkNvar[i][j].get(GRB_DoubleAttr_X);
                //     }
                // }

                /*********************************** follew reference *****************************************/

                // Update
                // update x0

                for (int j = 0; j < 3; j++)
                {
                    pos[j] = posNvar[0][j].get(GRB_DoubleAttr_X);
                    v[j] = vNvar[0][j].get(GRB_DoubleAttr_X);
                    a[j] = aNvar[0][j].get(GRB_DoubleAttr_X);
                    jerk[j] = jerkNvar[0][j].get(GRB_DoubleAttr_X);
                    pos_n[j] = posNvar[N - 1][j].get(GRB_DoubleAttr_X);
                }

                // whether close enough
                // if (check_update(pos, refposN[0]))
                if (check_update(pos, refposN[0], pos_n, refposN[N - 1]))
                {
                    ROS_INFO("---UPDATE---");
                    fail = 0;
                    iteration++;
                    // update reference
                    for (int i = 0; i < N; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            if (i + iteration < refNum)
                            {
                                refposN[i][j] = ref_nodes[i + iteration](j);
                            }
                            else
                            {
                                refposN[i][j] = ref_nodes[refNum - 1](j);
                            }
                        }
                    }
                    // update xN
                    for (int i = 0; i < N; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            posN[i][j] = posNvar[i][j].get(GRB_DoubleAttr_X);
                            vN[i][j] = vNvar[i][j].get(GRB_DoubleAttr_X);
                            aN[i][j] = aNvar[i][j].get(GRB_DoubleAttr_X);
                            jerkN[i][j] = jerkNvar[i][j].get(GRB_DoubleAttr_X);
                        }
                    }
                }
                else // update x0 = x1...
                {
                    fail++;
                    ROS_INFO("Agent %d Fails %d times", agent_id, fail);
                    if (fail >= N)
                    {
                        ROS_INFO("Agent %d : Solution Failure", agent_id);
                        for (int i = 0; i < 3; i++)
                        {
                            pos[i] = posN[N - 1][i];
                            v[i] = vN[N - 1][i];
                            a[i] = aN[N - 1][i];
                            jerk[i] = jerkN[N - 1][i];
                        }
                    }
                    else
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            pos[i] = posN[fail][i];
                            v[i] = vN[fail][i];
                            a[i] = aN[fail][i];
                            jerk[i] = jerkN[fail][i];
                        }
                    }
                }

                break;

            case GRB_INFEASIBLE:
                std::cout << "Model is infeasible." << std::endl;
                model.computeIIS();
                model.write("model_file.ilp");
                std::cout << "IIS:" << std::endl;
                for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); ++i)
                {
                    GRBConstr c = model.getConstr(i);
                    if (c.get(GRB_IntAttr_IISConstr) == 1)
                    {
                        std::cout << c.get(GRB_StringAttr_ConstrName) << std::endl;
                    }
                }

                fail++;
                ROS_INFO("Agent %d Fails %d times", agent_id, fail);
                if (fail >= N)
                {
                    ROS_INFO("Agent %d : Solution Failure", agent_id);
                    for (int i = 0; i < 3; i++)
                    {
                        pos[i] = posN[N - 1][i];
                        v[i] = vN[N - 1][i];
                        a[i] = aN[N - 1][i];
                        jerk[i] = jerkN[N - 1][i];
                    }
                }
                else
                {
                    for (int i = 0; i < 3; i++)
                    {
                        pos[i] = posN[fail][i];
                        v[i] = vN[fail][i];
                        a[i] = aN[fail][i];
                        jerk[i] = jerkN[fail][i];
                    }
                }

                break;
            case GRB_TIME_LIMIT:
                std::cout << "Time limit exceeded." << std::endl;

                // Update
                // update x0

                for (int j = 0; j < 3; j++)
                {
                    pos[j] = posNvar[0][j].get(GRB_DoubleAttr_X);
                    v[j] = vNvar[0][j].get(GRB_DoubleAttr_X);
                    a[j] = aNvar[0][j].get(GRB_DoubleAttr_X);
                    jerk[j] = jerkNvar[0][j].get(GRB_DoubleAttr_X);
                    pos_n[j] = posNvar[N - 1][j].get(GRB_DoubleAttr_X);
                }

                // whether close enough
                // if (check_update(pos, refposN[0]))
                if (check_update(pos, refposN[0], pos_n, refposN[N - 1]))
                {
                    ROS_INFO("---UPDATE---");
                    fail = 0;
                    iteration++;
                    // update reference
                    for (int i = 0; i < N; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            if (i + iteration < refNum)
                            {
                                refposN[i][j] = ref_nodes[i + iteration](j);
                            }
                            else
                            {
                                refposN[i][j] = ref_nodes[refNum - 1](j);
                            }
                        }
                    }
                    // update xN
                    for (int i = 0; i < N; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            posN[i][j] = posNvar[i][j].get(GRB_DoubleAttr_X);
                            vN[i][j] = vNvar[i][j].get(GRB_DoubleAttr_X);
                            aN[i][j] = aNvar[i][j].get(GRB_DoubleAttr_X);
                            jerkN[i][j] = jerkNvar[i][j].get(GRB_DoubleAttr_X);
                        }
                    }
                }
                else // update x0 = x1...
                {
                    fail++;
                    ROS_INFO("Agent %d Fails %d times", agent_id, fail);
                    if (fail >= N)
                    {
                        ROS_INFO("Agent %d : Solution Failure", agent_id);
                        for (int i = 0; i < 3; i++)
                        {
                            pos[i] = posN[N - 1][i];
                            v[i] = vN[N - 1][i];
                            a[i] = aN[N - 1][i];
                            jerk[i] = jerkN[N - 1][i];
                        }
                    }
                    else
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            pos[i] = posN[fail][i];
                            v[i] = vN[fail][i];
                            a[i] = aN[fail][i];
                            jerk[i] = jerkN[fail][i];
                        }
                    }
                }

                break;
            case GRB_UNBOUNDED:
                std::cout << "Model is unbounded." << std::endl;
                break;
            case GRB_ERROR_NUMERIC:
                std::cout << "Numeric issues occurred." << std::endl;
                break;
            case GRB_INTERRUPTED:
                std::cout << "Optimization was interrupted." << std::endl;
                break;

            case GRB_ITERATION_LIMIT:
                std::cout << "Iteration limit exceeded." << std::endl;
                break;
            case GRB_NODE_LIMIT:
                std::cout << "Node limit exceeded." << std::endl;
                break;
            case GRB_USER_OBJ_LIMIT:
                std::cout << "User-specified objective limit reached." << std::endl;
                break;
            case GRB_BATCH_STATUS_UNKNOWN:
                std::cout << "Unknown status." << std::endl;
                break;
            default:
                std::cout << "Unexpected status code: " << status << std::endl;
            }

            // // print value of variables
            std::cout << "AGENT_ID = " << agent_id << "  NUM = " << num << "  ITERATION = " << iteration << std::endl;
            std::cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl
                      << std::endl;
            //           << std::endl
            //           << std::endl;

            // Record

            // 创建一个随机数引擎
            std::random_device rd;
            std::default_random_engine generator(rd());
            // 创建正态分布对象
            std::normal_distribution<double> distribution_pos(0, qp);
            std::normal_distribution<double> distribution_v(0, qv);
            double nose;

            // Eigen::Vector3d coord;
            // coord(0) = pos[0] + distribution_pos(generator);
            // coord(1) = pos[1] + distribution_pos(generator);
            // coord(2) = pos[2] + distribution_pos(generator);
            // pos_nodes.push_back(coord);

            // Eigen::Vector3d speed;
            // speed(0) = v[0] + distribution_v(generator);
            // speed(1) = v[1] + distribution_v(generator);
            // speed(2) = v[2] + distribution_v(generator);
            // v_nodes.push_back(speed);

            Eigen::Vector3d coord;
            coord(0) = pos[0];
            coord(1) = pos[1];
            coord(2) = pos[2];
            pos_nodes.push_back(coord);

            Eigen::Vector3d speed;
            speed(0) = v[0];
            speed(1) = v[1];
            speed(2) = v[2];
            v_nodes.push_back(speed);

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    posN_nodes[j](i) = posN[j][i];
                }
            }

            num++;
        }
        catch (GRBException e)
        {
            std::cout << "Error code " << e.getErrorCode() << " : " << e.getMessage() << std::endl;
        }
        catch (...)
        {
            std::cout << "Exception during optimization" << std::endl;
        }
        // ROS_INFO("NUM = %d, ITERATION = %d", num, iteration);
        // ROS_INFO("x = %f, y = %f, z = %f", pos[0], pos[1], pos[2]);
        // ROS_INFO("vx = %f, vy = %f, vz = %f", v[0], v[1], v[2]);
        // ROS_INFO("ax = %f, ay = %f, az = %f", a[0], a[1], a[2]);
        // ROS_INFO("jx = %f, jy = %f, jz = %f", jerk[0], jerk[1], jerk[2]);
        // ROS_INFO("refx = %f, refy = %f, refz = %f", refposN[0][0], refposN[0][1], refposN[0][2]);

        // std::cout << "agent " << agent_id << " posN: " << std::endl;
        // for (int i = 0; i < 3; i++)
        // {
        //     for (int j = 0; j < N; j++)
        //     {
        //         std::cout << posN[j][i] << "  ";
        //     }
        //     std::cout << std::endl;
        // }

        // std::cout << std::endl;
        // std::cout << "refN:" << std::endl;
        // for (int i = 0; i < 3; i++)
        // {
        //     for (int j = 0; j < N; j++)
        //     {
        //         std::cout << refposN[j][i] << "  ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << std::endl;
    }
    else
    {
        ROS_INFO("Agent %d Planning Failure !", agent_id);
        flag_failure = true;
    }
}
