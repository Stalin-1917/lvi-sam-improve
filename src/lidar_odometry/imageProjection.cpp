#include "utility.h"
#include "lvi_sam/cloud_info.h"
/**
 * Velodyne点云结构，变量名XYZIRT是每个变量的首字母
*/
struct VelodynePointXYZIRT
{
    PCL_ADD_POINT4D
    PCL_ADD_INTENSITY;
    uint16_t ring;
    float time;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
} EIGEN_ALIGN16;
// 注册为PCL点云格式
POINT_CLOUD_REGISTER_POINT_STRUCT(VelodynePointXYZIRT,
                                  (float, x, x)(float, y, y)(float, z, z)(float, intensity, intensity)(uint16_t, ring, ring)(float, time, time))

struct OusterPointXYZIRT
{
    PCL_ADD_POINT4D;
    float intensity;
    uint32_t t;
    uint16_t reflectivity;
    uint8_t ring;
    uint16_t noise;
    uint32_t range;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
} EIGEN_ALIGN16;
POINT_CLOUD_REGISTER_POINT_STRUCT(OusterPointXYZIRT,
                                  (float, x, x)(float, y, y)(float, z, z)(float, intensity, intensity)(uint32_t, t, t)(uint16_t, reflectivity, reflectivity)(uint8_t, ring, ring)(uint16_t, noise, noise)(uint32_t, range, range))

// Use the Velodyne point format as a common representation
using PointXYZIRT = VelodynePointXYZIRT;

const int queueLength = 2000;
class AKFmanager
{
    private:
        int type;
        int state;
        Eigen::Matrix<double, 12, 1> x;
        Eigen::Matrix<double, 12, 1> z;
        Eigen::Matrix<double, 12,12> H;
        Eigen::Matrix<double, 12,12> A;
        Eigen::Matrix<double, 12,12> K;
        Eigen::Matrix<double, 12,12> P;
        Eigen::Matrix<double, 12,12> P_1;
        Eigen::Matrix<double, 12,12> Q;
        Eigen::Matrix<double, 12,12> R;
        Eigen::Matrix<double, 12, 1> epsilon;
        Eigen::Matrix<double, 12 ,1> delta;
        float factor; 
    public:
        
        AKFmanager(Eigen::Matrix<double,12,12> A_init,Eigen::Matrix<double,12,12> H_init,Eigen::Matrix<double, 12,12> P_init,Eigen::Matrix<double, 12,12> Q_init,Eigen::Matrix<double, 12,12> R_init,float alpha)
        {
            state = 0;
            A = A_init;
            H = H_init;
            P = P_init;
            Q = Q_init;
            R = R_init;
            
            factor = alpha;
        }
        AKFmanager(float alpha)
        {
            ROS_INFO("AKF_START");
            factor = alpha;
            state = 0;
            A = Eigen::MatrixXd::Identity(12,12);
            H = Eigen::MatrixXd::Identity(12,12);
            P = 1000000 * Eigen::MatrixXd::Identity(12,12);
            Q = 0.01 * Eigen::MatrixXd::Identity(12,12);
            R = 0.01 * Eigen::MatrixXd::Identity(12,12);
            
        }
        AKFmanager()
        {
            ROS_INFO("AKF_START");
            factor = 0.9;
            state = 0;
            A = Eigen::MatrixXd::Identity(12,12);
            H = Eigen::MatrixXd::Identity(12,12);
            P = 1000000 * Eigen::MatrixXd::Identity(12,12);
            Q = 0.01 * Eigen::MatrixXd::Identity(12,12);
            R = 0.01 * Eigen::MatrixXd::Identity(12,12);
            
            
        }
        ~AKFmanager() {}
        void input(vector<float> sensor)
        {
            for (int i = 0 ; i < sensor.size();i++)
            {
                if (i >= 12)
                    return;
                z[i] = sensor.at(i);
            }
        }
        void input(Eigen::Matrix<double,12,1> sensor)
        {
            //ROS_INFO("input!");
            z = sensor;
        }
        void input(vector<float> odom,vector<float> odomIncre)
        {
            for (int i = 0 ; i < odom.size();i++)
            {
                if (i >= 6)
                    break;
                z[i] = odom.at(i);
            }
            for (int i = 0 ; i < odomIncre.size();i++)
            {
                if (i >= 6)
                    break;
                z[6 + i] = odomIncre.at(i);
            }
        }
        
        Eigen::Matrix<double, 12, 1> forward()
        {
            //ROS_INFO("step_0");
            if (state == 0)
            {
                //x = H.inverse()*z;
                x = z;
                state = 1;
                return x;
            }
           // ROS_INFO("step_1");
            x = A * x;
            //ROS_INFO("step_2");
            delta = z - H * x;
            //ROS_INFO("step_3");
            P = A * P * A.transpose() + Q;
            //ROS_INFO("step_4");
            K = P * H.transpose()*(H * P * H.transpose() + R).inverse();
            //ROS_INFO("step_5");
            x = x + K * delta;
            //ROS_INFO("step_6");
            epsilon = z - H * x;
            //ROS_INFO("step_7");

            R = factor * R + (1-factor)*(epsilon*epsilon.transpose() + H * P * H.transpose());
            P = (Eigen::MatrixXd::Identity(12,12) - K * H) * P;
            //ROS_INFO("step_8");
            Q = factor * Q + (1-factor)*(K * delta * delta.transpose() * K.transpose());
            
            return x;
        }
        void update_A(float dtime)
        {
            A.block(0,6,6,6) = dtime * Eigen::MatrixXd::Identity(6,6);
        }
        Eigen::Matrix<double, 12, 1> getResult(void)
        {
            return x;
        }
        void reset_system()
        {
            state = 0;
            Q = 0.01 * Eigen::MatrixXd::Identity(12,12);
            R = 0.01 * Eigen::MatrixXd::Identity(12,12);

        }
};
class ImageProjection : public ParamServer 
{
private:
    AKFmanager AKF_imu_vin;
    std::mutex imuLock;
    std::mutex vinsOdomLock;
    std::mutex imuOdomLock;

    ros::Subscriber subLaserCloud;
    ros::Publisher pubLaserCloud;

    ros::Publisher pubExtractedCloud;
    ros::Publisher pubLaserCloudInfo;

    ros::Subscriber subImu;
    std::deque<sensor_msgs::Imu> imuQueue;

    ros::Subscriber subVinsOdom;
    ros::Subscriber subImuOdom;
    std::deque<nav_msgs::Odometry> vinsOdomQueue;
    std::deque<nav_msgs::Odometry> imuOdomQueue;

    std::deque<sensor_msgs::PointCloud2> cloudQueue;
    sensor_msgs::PointCloud2 currentCloudMsg;

    double *imuTime = new double[queueLength];
    double *imuRotX = new double[queueLength];
    double *imuRotY = new double[queueLength];
    double *imuRotZ = new double[queueLength];

    int imuPointerCur;
    bool firstPointFlag;
    Eigen::Affine3f transStartInverse;

    pcl::PointCloud<PointXYZIRT>::Ptr laserCloudIn;
    pcl::PointCloud<OusterPointXYZIRT>::Ptr tmpOusterCloudIn;
    pcl::PointCloud<PointType>::Ptr fullCloud;
    pcl::PointCloud<PointType>::Ptr extractedCloud;
    
    int deskewFlag;
    cv::Mat rangeMat;

    bool odomDeskewFlag;
    bool vinsDeskewFlag;
    bool AKF_deskewFlag;
    bool imuDeskewFlag;
    float odomIncreX;
    float odomIncreY;
    float odomIncreZ;

    float odomIncreImu[6];
    float odomIncreVin[6];
    float odomImu[6];
    float odomVin[6];
    float odomIncre[6];

    nav_msgs::Odometry startImuOdomMsg;
    nav_msgs::Odometry endImuOdomMsg;

    nav_msgs::Odometry startVinOdomMsg;
    nav_msgs::Odometry endVinOdomMsg;

    lvi_sam::cloud_info cloudInfo;
    double timeScanCur;
    double timeScanEnd;
    std_msgs::Header cloudHeader;

    vector<int> columnIdnCountVec;
    


public:
    ImageProjection() : deskewFlag(0)
    {
        // 订阅原始imu数据
        subImu = nh.subscribe<sensor_msgs::Imu>(imuTopic, 2000, &ImageProjection::imuHandler, this, ros::TransportHints().tcpNoDelay());
        //! 重要：VIO发来的里程计消息，会作为后端点云配准的位姿初值
        subVinsOdom = nh.subscribe<nav_msgs::Odometry>(PROJECT_NAME + "/vins/odometry/imu_propagate_ros", 2000, &ImageProjection::vinsOdometryHandler, this, ros::TransportHints().tcpNoDelay());
        // 订阅imu里程计，由imuPreintegration积分计算得到的每时刻imu位姿
    
        subImuOdom = nh.subscribe<nav_msgs::Odometry>(odomTopic + "_incremental", 2000, &ImageProjection::imuOdometryHandler, this, ros::TransportHints().tcpNoDelay());
        // 订阅原始lidar数据
    
        subLaserCloud = nh.subscribe<sensor_msgs::PointCloud2>(pointCloudTopic, 5, &ImageProjection::cloudHandler, this, ros::TransportHints().tcpNoDelay());
        // 发布当前激光帧运动畸变校正后的点云，有效点
        pubExtractedCloud = nh.advertise<sensor_msgs::PointCloud2>(PROJECT_NAME + "/lidar/deskew/cloud_deskewed", 5);
        // 发布当前激光帧运动畸变校正后的点云信息
        pubLaserCloudInfo = nh.advertise<lvi_sam::cloud_info>(PROJECT_NAME + "/lidar/deskew/cloud_info", 5);

        allocateMemory();
        resetParameters();

        pcl::console::setVerbosityLevel(pcl::console::L_ERROR);
    }

    void allocateMemory()
    {
        laserCloudIn.reset(new pcl::PointCloud<PointXYZIRT>());
        tmpOusterCloudIn.reset(new pcl::PointCloud<OusterPointXYZIRT>());
        fullCloud.reset(new pcl::PointCloud<PointType>());
        extractedCloud.reset(new pcl::PointCloud<PointType>());

        fullCloud->points.resize(N_SCAN * Horizon_SCAN);

        cloudInfo.startRingIndex.assign(N_SCAN, 0);
        cloudInfo.endRingIndex.assign(N_SCAN, 0);

        cloudInfo.pointColInd.assign(N_SCAN * Horizon_SCAN, 0);
        cloudInfo.pointRange.assign(N_SCAN * Horizon_SCAN, 0);
        for(int i = 0;i<6;i++)
        {
            odomIncreImu[i] = 0;
            odomIncreVin[i] = 0;
            odomImu[i] = 0;
            odomVin[i] = 0;
            odomIncre[i] = 0;
        }

        resetParameters();
    }

    void resetParameters()
    {
        laserCloudIn->clear();
        extractedCloud->clear();
        // reset range matrix for range image projection
        rangeMat = cv::Mat(N_SCAN, Horizon_SCAN, CV_32F, cv::Scalar::all(FLT_MAX));

        imuPointerCur = 0;
        firstPointFlag = true;
        odomDeskewFlag = false;

        AKF_deskewFlag = false;
        for (int i = 0; i < queueLength; ++i)
        {
            imuTime[i] = 0;
            imuRotX[i] = 0;
            imuRotY[i] = 0;
            imuRotZ[i] = 0;
        }

        columnIdnCountVec.assign(N_SCAN, 0);
    }

    ~ImageProjection() {}

    void imuHandler(const sensor_msgs::Imu::ConstPtr &imuMsg)
    {
        sensor_msgs::Imu thisImu = imuConverter(*imuMsg);
        std::lock_guard<std::mutex> lock1(imuLock);
        imuQueue.push_back(thisImu);

        // debug IMU data
        // cout << std::setprecision(6);
        // cout << "IMU acc: " << endl;
        // cout << "x: " << thisImu.linear_acceleration.x << 
        //       ", y: " << thisImu.linear_acceleration.y << 
        //       ", z: " << thisImu.linear_acceleration.z << endl;
        // cout << "IMU gyro: " << endl;
        // cout << "x: " << thisImu.angular_velocity.x << 
        //       ", y: " << thisImu.angular_velocity.y << 
        //       ", z: " << thisImu.angular_velocity.z << endl;
        // double imuRoll, imuPitch, imuYaw;
        // tf::Quaternion orientation;
        // tf::quaternionMsgToTF(thisImu.orientation, orientation);
        // tf::Matrix3x3(orientation).getRPY(imuRoll, imuPitch, imuYaw);
        // cout << "IMU roll pitch yaw: " << endl;
        // cout << "roll: " << imuRoll << ", pitch: " << imuPitch << ", yaw: " << imuYaw << endl << endl;

    }
/**
 * 订阅视觉里程计
*/
    void vinsOdometryHandler(const nav_msgs::Odometry::ConstPtr &odometryMsg)
    {
        std::lock_guard<std::mutex> lock2(vinsOdomLock);
        vinsOdomQueue.push_back(*odometryMsg);
    }
/**
 * 订阅imu里程计，由imuPreintegration积分计算得到的每时刻imu位姿
*/
    void imuOdometryHandler(const nav_msgs::Odometry::ConstPtr &odometryMsg)
    {
        std::lock_guard<std::mutex> lock2(imuOdomLock);
        imuOdomQueue.push_back(*odometryMsg);
    }

    void cloudHandler(const sensor_msgs::PointCloud2ConstPtr &laserCloudMsg)
    {
        if (!cachePointCloud(laserCloudMsg))
            return;

        if (!deskewInfo())
            return;

        projectPointCloud();

        cloudExtraction();

        publishClouds();

        resetParameters();
    }
/*添加一帧激光点云到队列，取出最早一帧作为当前帧，计算起止时间戳，检查数据有效性*/
    bool cachePointCloud(const sensor_msgs::PointCloud2ConstPtr &laserCloudMsg)
    {
        // cache point cloud
        cloudQueue.push_back(*laserCloudMsg);
        if (cloudQueue.size() <= 2)
            return false;

        // convert cloud
        currentCloudMsg = std::move(cloudQueue.front());
        cloudQueue.pop_front();
        if (sensor == SensorType::VELODYNE || sensor == SensorType::LIVOX)
        {
            pcl::moveFromROSMsg(currentCloudMsg, *laserCloudIn);
        }
        else if (sensor == SensorType::OUSTER)
        {
            // Convert to Velodyne format
            pcl::moveFromROSMsg(currentCloudMsg, *tmpOusterCloudIn);
            laserCloudIn->points.resize(tmpOusterCloudIn->size());
            laserCloudIn->is_dense = tmpOusterCloudIn->is_dense;
            for (size_t i = 0; i < tmpOusterCloudIn->size(); i++)
            {
                auto &src = tmpOusterCloudIn->points[i];
                auto &dst = laserCloudIn->points[i];
                dst.x = src.x;
                dst.y = src.y;
                dst.z = src.z;
                dst.intensity = src.intensity;
                dst.ring = src.ring;
                dst.time = src.t * 1e-9f;
            }
        }
        else
        {
            ROS_ERROR_STREAM("Unknown sensor type: " << int(sensor));
            ros::shutdown();
        }

        // get timestamp
        cloudHeader = currentCloudMsg.header;
        timeScanCur = cloudHeader.stamp.toSec();
        timeScanEnd = timeScanCur + laserCloudIn->points.back().time;

        // check dense flag
        if (laserCloudIn->is_dense == false)
        {
            ROS_ERROR("Point cloud is not in dense format, please remove NaN points first!");
            ros::shutdown();
        }

        // check ring channel
        static int ringFlag = 0;
        if (ringFlag == 0)
        {
            ringFlag = -1;
            for (int i = 0; i < (int)currentCloudMsg.fields.size(); ++i)
            {
                if (currentCloudMsg.fields[i].name == "ring")
                {
                    ringFlag = 1;
                    break;
                }
            }
            if (ringFlag == -1)
            {
                ROS_ERROR("Point cloud ring channel not available, please configure your point cloud data!");
                ros::shutdown();
            }
        }

        // check point time
        if (deskewFlag == 0)
        {
            deskewFlag = -1;
            for (auto &field : currentCloudMsg.fields)
            {
                if (field.name == "time" || field.name == "t")
                {
                    deskewFlag = 1;
                    break;
                }
            }
            if (deskewFlag == -1)
                ROS_WARN("Point cloud timestamp not available, deskew function disabled, system will drift significantly!");
        }

        return true;
    }
/**
 * 当前帧起止时刻对应的imu数据、imu里程计数据处理
*/
    bool deskewInfo()
    {
        std::lock_guard<std::mutex> lock1(imuLock);

        // make sure IMU data available for the scan
        if (imuQueue.empty() || imuQueue.front().header.stamp.toSec() > timeScanCur || imuQueue.back().header.stamp.toSec() < timeScanEnd)
        {
            ROS_DEBUG("Waiting for IMU data ...");
            return false;
        }

        imuDeskewInfo();

        {
            std::lock_guard<std::mutex> lock2(vinsOdomLock);
            vinsOdomDeskewInfo();
            //vinsOdomDeskewInfo_old();
        }

        {
            std::lock_guard<std::mutex> lock2(imuOdomLock);
            imuOdomDeskewInfo();
            //imuOdomDeskewInfo_old();
        }
        if (!odomDeskewFlag)
        {
            cloudInfo.initialGuessX = 0;
            cloudInfo.initialGuessY = 0;
            cloudInfo.initialGuessZ = 0;
            cloudInfo.initialGuessRoll = 0;
            cloudInfo.initialGuessPitch = 0;
            cloudInfo.initialGuessYaw = 0;
            odomIncreX = 0;
            odomIncreY = 0;
            odomIncreZ = 0;
            return true;
        }
            
        // changed codes
        //使用imu得到的增量和vins得到的绝对坐标

        if (cloudInfo.vinsOdomAvailable == false)
        {

            double roll,pitch,yaw;
            tf::Quaternion orientation;
            tf::quaternionMsgToTF(startImuOdomMsg.pose.pose.orientation, orientation);
            tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);
            // Initial guess used in mapOptimization
            //; 这里就使用imu odom作为后端优化的初值
            cloudInfo.initialGuessX = odomImu[0];
            cloudInfo.initialGuessY = odomImu[1];
            cloudInfo.initialGuessZ = odomImu[2];
            cloudInfo.initialGuessRoll = odomImu[3];
            cloudInfo.initialGuessPitch = odomImu[4];
            cloudInfo.initialGuessYaw = odomImu[5];
            for (int i = 0 ; i < 6 ; i ++)
            {
                odomIncre[i] = odomIncreImu[i];
            }
        }
        else
        {
            
            AKF_deskewFlag = false;
            Eigen::Matrix<double,12,1> sensor;
            Eigen::Matrix<double,12,1> AKF_results;
            tf::Quaternion orientation;
            double roll,pitch,yaw;
            tf::quaternionMsgToTF(startImuOdomMsg.pose.pose.orientation, orientation);
            tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);

            for (int i = 0 ; i < 6 ; i ++)
            {
                sensor[i] = odomVin[i];
                sensor[6 + i] = odomIncreImu[i];
            }
            AKF_imu_vin.update_A(laserCloudIn->points.back().time);
            AKF_imu_vin.input(sensor);
            AKF_results = AKF_imu_vin.forward();

            Eigen::Matrix<double,12,1> Test_results;
            if (vinsDeskewFlag == true)
            {
                for (int i = 0;i < 6;i++)
                {
                    Test_results[i] = odomVin[i];
                    Test_results[i + 6] = odomIncreVin[i];
                }
                if ((Test_results.block(0,0,6,1) - AKF_results.block(0,0,6,1)).norm () < 0.06)//源程序运行误差不大，滤波结果偏差太大肯定是失败的
                {
                    cloudInfo.initialGuessX = AKF_results[0];
                    cloudInfo.initialGuessY = AKF_results[1];
                    cloudInfo.initialGuessZ = AKF_results[2];
                    cloudInfo.initialGuessRoll = AKF_results[3];
                    cloudInfo.initialGuessPitch = AKF_results[4];
                    cloudInfo.initialGuessYaw = AKF_results[5];
                    for (int i = 0;i < 6;i++)
                    {
                        odomIncre[i] = AKF_results[6 + i];
                    }
                    ROS_INFO("USE AKF RESULTS");
                }
                    

                else
                {
                    ROS_INFO("USE VIN RESULTs");
                    AKF_imu_vin.reset_system();
                    if (cloudInfo.vinsOdomAvailable == true)
                    {
                        cloudInfo.initialGuessX = odomVin[0];
                        cloudInfo.initialGuessY = odomVin[1];
                        cloudInfo.initialGuessZ = odomVin[2];
                        cloudInfo.initialGuessRoll = odomVin[3];
                        cloudInfo.initialGuessPitch = odomVin[4];
                        cloudInfo.initialGuessYaw = odomVin[5];
                    }
                    if (vinsDeskewFlag)
                    {
                        for (int i = 0 ; i < 6 ; i ++)
                        {
                            odomIncre[i] = odomIncreVin[i];
                        }
                    }
                    }
            }
            else
            {
                ROS_INFO("USE IMU RESULTs");
                if (cloudInfo.imuOdomAvailable == true)
                {
                    cloudInfo.initialGuessX = odomImu[0];
                    cloudInfo.initialGuessY = odomImu[1];
                    cloudInfo.initialGuessZ = odomImu[2];
                    cloudInfo.initialGuessRoll = odomImu[3];
                    cloudInfo.initialGuessPitch = odomImu[4];
                    cloudInfo.initialGuessYaw = odomImu[5];
                
                }
                if (  imuDeskewFlag)
                {
                    for (int i = 0 ; i < 6 ; i ++)
                    {
                        odomIncre[i] = odomIncreImu[i];
                    }
                }
                AKF_imu_vin.reset_system();
                
            }
            AKF_deskewFlag = true;
        }
        odomIncreX = odomIncre[0];
        odomIncreY = odomIncre[1];
        odomIncreZ = odomIncre[2];
        
        
        //AKF output cloudInfo.initialGuess and odomIncreX, odomIncreY, odomIncreZ
        return true;
    }
    void imuDeskewInfo()
    {
        cloudInfo.imuAvailable = false;

        while (!imuQueue.empty())
        {
            if (imuQueue.front().header.stamp.toSec() < timeScanCur - 0.01)
                imuQueue.pop_front();
            else
                break;
        }

        if (imuQueue.empty())
            return;

        imuPointerCur = 0;

        for (int i = 0; i < (int)imuQueue.size(); ++i)
        {
            sensor_msgs::Imu thisImuMsg = imuQueue[i];
            double currentImuTime = thisImuMsg.header.stamp.toSec();

            // get roll, pitch, and yaw estimation for this scan
            if (currentImuTime <= timeScanCur)
                imuRPY2rosRPY(&thisImuMsg, &cloudInfo.imuRollInit, &cloudInfo.imuPitchInit, &cloudInfo.imuYawInit);

            if (currentImuTime > timeScanEnd + 0.01)
                break;

            if (imuPointerCur == 0)
            {
                imuRotX[0] = 0;
                imuRotY[0] = 0;
                imuRotZ[0] = 0;
                imuTime[0] = currentImuTime;
                ++imuPointerCur;
                continue;
            }

            // get angular velocity
            double angular_x, angular_y, angular_z;
            imuAngular2rosAngular(&thisImuMsg, &angular_x, &angular_y, &angular_z);

            // integrate rotation
            double timeDiff = currentImuTime - imuTime[imuPointerCur - 1];
            imuRotX[imuPointerCur] = imuRotX[imuPointerCur - 1] + angular_x * timeDiff;
            imuRotY[imuPointerCur] = imuRotY[imuPointerCur - 1] + angular_y * timeDiff;
            imuRotZ[imuPointerCur] = imuRotZ[imuPointerCur - 1] + angular_z * timeDiff;
            imuTime[imuPointerCur] = currentImuTime;
            ++imuPointerCur;
        }

        --imuPointerCur;

        if (imuPointerCur <= 0)
            return;

        cloudInfo.imuAvailable = true;
    }

    void vinsOdomDeskewInfo_old()
    {
        cloudInfo.vinsOdomAvailable = false;

        while (!vinsOdomQueue.empty())
        {
            if (vinsOdomQueue.front().header.stamp.toSec() < timeScanCur - 0.01)
                vinsOdomQueue.pop_front();
            else
                break;
        }

        if (vinsOdomQueue.empty())
            return;

        if (vinsOdomQueue.front().header.stamp.toSec() > timeScanCur)
            return;

        // get start odometry at the beinning of the scan
        nav_msgs::Odometry startOdomMsg;

        for (int i = 0; i < (int)vinsOdomQueue.size(); ++i)
        {
            startOdomMsg = vinsOdomQueue[i];

            if (ROS_TIME(&startOdomMsg) < timeScanCur)
                continue;
            else
                break;
        }

        tf::Quaternion orientation;
        tf::quaternionMsgToTF(startOdomMsg.pose.pose.orientation, orientation);

        double roll, pitch, yaw;
        tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);

        // Initial guess used in mapOptimization
        //! 重要：这里会把前端vins发来的位姿作为后端scan-to-map匹配的初值，所以vins发来的必须是T_world_lidar
        cloudInfo.initialGuessX = startOdomMsg.pose.pose.position.x;
        cloudInfo.initialGuessY = startOdomMsg.pose.pose.position.y;
        cloudInfo.initialGuessZ = startOdomMsg.pose.pose.position.z;
        cloudInfo.initialGuessRoll = roll;
        cloudInfo.initialGuessPitch = pitch;
        cloudInfo.initialGuessYaw = yaw;
        //; vins里程计重启id，在计算后端优化的初值时会使用
        cloudInfo.vinsOdomResetId = (int)round(startOdomMsg.pose.covariance[0]);
        cloudInfo.vinsOdomAvailable = true;

        // get end odometry at the end of the scan
        odomDeskewFlag = false;

        if (vinsOdomQueue.back().header.stamp.toSec() < timeScanEnd)
            return;

        nav_msgs::Odometry endOdomMsg;

        for (int i = 0; i < (int)vinsOdomQueue.size(); ++i)
        {
            endOdomMsg = vinsOdomQueue[i];

            if (ROS_TIME(&endOdomMsg) < timeScanEnd)
                continue;
            else
                break;
        }

        //; 要保证前后的vins odom数据的id是一样的，即vins odom没有经过复位操作
        if (int(round(startOdomMsg.pose.covariance[0])) != int(round(endOdomMsg.pose.covariance[0])))
            return;

        Eigen::Affine3f transBegin = pcl::getTransformation(startOdomMsg.pose.pose.position.x, startOdomMsg.pose.pose.position.y, startOdomMsg.pose.pose.position.z, roll, pitch, yaw);

        tf::quaternionMsgToTF(endOdomMsg.pose.pose.orientation, orientation);
        tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);
        Eigen::Affine3f transEnd = pcl::getTransformation(endOdomMsg.pose.pose.position.x, endOdomMsg.pose.pose.position.y, endOdomMsg.pose.pose.position.z, roll, pitch, yaw);

        Eigen::Affine3f transBt = transBegin.inverse() * transEnd;

        float rollIncre, pitchIncre, yawIncre;
        pcl::getTranslationAndEulerAngles(transBt, odomIncreX, odomIncreY, odomIncreZ, rollIncre, pitchIncre, yawIncre);

        odomDeskewFlag = true;
    }
//这里可以改成融合滤波

    void imuOdomDeskewInfo_old()

    {
        cloudInfo.imuOdomAvailable = false;
        nav_msgs::Odometry startOdomMsg;
        tf::Quaternion orientation;
        double roll, pitch, yaw;
        
        while (!imuOdomQueue.empty())
        {
            if (imuOdomQueue.front().header.stamp.toSec() < timeScanCur - 0.01)
                imuOdomQueue.pop_front();
            else
                break;
        }

        if (imuOdomQueue.empty())
            return;

        if (imuOdomQueue.front().header.stamp.toSec() > timeScanCur)
            return;

        //? add: 如果vins odom不可用，则寻找imu odom
        if (cloudInfo.vinsOdomAvailable == false)
        {
            // get start odometry at the beinning of the scan

            for (int i = 0; i < (int)imuOdomQueue.size(); ++i)
            {
                startOdomMsg = imuOdomQueue[i];

                if (ROS_TIME(&startOdomMsg) < timeScanCur)
                    continue;
                else
                    break;
            }

            tf::quaternionMsgToTF(startOdomMsg.pose.pose.orientation, orientation);
            tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);

            // Initial guess used in mapOptimization
            //; 这里就使用imu odom作为后端优化的初值
            cloudInfo.initialGuessX = startOdomMsg.pose.pose.position.x;
            cloudInfo.initialGuessY = startOdomMsg.pose.pose.position.y;
            cloudInfo.initialGuessZ = startOdomMsg.pose.pose.position.z;
            cloudInfo.initialGuessRoll = roll;
            cloudInfo.initialGuessPitch = pitch;
            cloudInfo.initialGuessYaw = yaw;
            //; imu里程计重启id，在计算后端优化的初值时会使用
            cloudInfo.imuOdomResetId = (int)round(startOdomMsg.pose.covariance[0]);
            cloudInfo.imuOdomAvailable = true;
        }

        //? add: 同理，如果vins odom的增量平移变换不可用，则寻找imu odom的增量平移变换        
        if (odomDeskewFlag == false)
        {
            // get end odometry at the end of the scan
            if (imuOdomQueue.back().header.stamp.toSec() < timeScanEnd)
                return;

            nav_msgs::Odometry endOdomMsg;

            for (int i = 0; i < (int)imuOdomQueue.size(); ++i)
            {
                endOdomMsg = imuOdomQueue[i];

                if (ROS_TIME(&endOdomMsg) < timeScanEnd)
                    continue;
                else
                    break;
            }

            //; 要保证前后的imu odom数据的id是一样的，即imu odom没有经过复位操作
            if (int(round(startOdomMsg.pose.covariance[0])) != int(round(endOdomMsg.pose.covariance[0])))
                return;

            Eigen::Affine3f transBegin = pcl::getTransformation(startOdomMsg.pose.pose.position.x, startOdomMsg.pose.pose.position.y, startOdomMsg.pose.pose.position.z, roll, pitch, yaw);
            tf::quaternionMsgToTF(endOdomMsg.pose.pose.orientation, orientation);
            double roll, pitch, yaw;
            tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);
            Eigen::Affine3f transEnd = pcl::getTransformation(endOdomMsg.pose.pose.position.x, endOdomMsg.pose.pose.position.y, endOdomMsg.pose.pose.position.z, roll, pitch, yaw);

            Eigen::Affine3f transBt = transBegin.inverse() * transEnd;

            float rollIncre, pitchIncre, yawIncre;
            pcl::getTranslationAndEulerAngles(transBt, odomIncreX, odomIncreY, odomIncreZ, rollIncre, pitchIncre, yawIncre);

            odomDeskewFlag = true;
        }
    }
    void vinsOdomDeskewInfo()
    {
        cloudInfo.vinsOdomAvailable = false;
        double roll, pitch, yaw;
        while (!vinsOdomQueue.empty())
        {
            if (vinsOdomQueue.front().header.stamp.toSec() < timeScanCur - 0.01)
                vinsOdomQueue.pop_front();
            else
                break;
        }

        if (vinsOdomQueue.empty())
        {
            ROS_INFO("empty Vins");
            return;
        }
        

        if (vinsOdomQueue.front().header.stamp.toSec() > timeScanCur)
        {
            ROS_INFO("vins_late");
            return;
        }
            

        // get start odometry at the beinning of the scan
        nav_msgs::Odometry startOdomMsg;

        for (int i = 0; i < (int)vinsOdomQueue.size(); ++i)
        {
            startVinOdomMsg = vinsOdomQueue[i];

            if (ROS_TIME(&startVinOdomMsg) < timeScanCur)
                continue;
            else
                break;
        }

        tf::Quaternion orientation;
        tf::quaternionMsgToTF(startVinOdomMsg.pose.pose.orientation, orientation);

        tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);
        odomVin[0] = startVinOdomMsg.pose.pose.position.x;
        odomVin[1] = startVinOdomMsg.pose.pose.position.y;
        odomVin[2] = startVinOdomMsg.pose.pose.position.z;
        odomVin[3] = roll;
        odomVin[4] = pitch;
        odomVin[5] = yaw;
        cloudInfo.vinsOdomAvailable = true;
        // get end odometry at the end of the scan
        vinsDeskewFlag = false;
        odomDeskewFlag = false;
        if (vinsOdomQueue.back().header.stamp.toSec() < timeScanEnd)
            return;

        nav_msgs::Odometry endOdomMsg;

        for (int i = 0; i < (int)vinsOdomQueue.size(); ++i)
        {
            endVinOdomMsg = vinsOdomQueue[i];

            if (ROS_TIME(&endVinOdomMsg) < timeScanEnd)
                continue;
            else
                break;
        }

        //; 要保证前后的vins odom数据的id是一样的，即vins odom没有经过复位操作
        if (int(round(startVinOdomMsg.pose.covariance[0])) != int(round(endVinOdomMsg.pose.covariance[0])))
            return;

        Eigen::Affine3f transBegin = pcl::getTransformation(startVinOdomMsg.pose.pose.position.x, startVinOdomMsg.pose.pose.position.y, startVinOdomMsg.pose.pose.position.z, roll, pitch, yaw);

        tf::quaternionMsgToTF(endVinOdomMsg.pose.pose.orientation, orientation);
        tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);
        Eigen::Affine3f transEnd = pcl::getTransformation(endVinOdomMsg.pose.pose.position.x, endVinOdomMsg.pose.pose.position.y, endVinOdomMsg.pose.pose.position.z, roll, pitch, yaw);

        Eigen::Affine3f transBt = transBegin.inverse() * transEnd;

        //float rollIncre, pitchIncre, yawIncre;

        pcl::getTranslationAndEulerAngles(transBt, odomIncreVin[0], odomIncreVin[1], odomIncreVin[2], odomIncreVin[3], odomIncreVin[4], odomIncreVin[5]);
        vinsDeskewFlag = true;
        odomDeskewFlag = true;
    }
    
    void imuOdomDeskewInfo()
    {

        cloudInfo.imuAvailable = false;
        nav_msgs::Odometry startOdomMsg;
        tf::Quaternion orientation;
        double roll, pitch, yaw;

        while (!imuOdomQueue.empty())
        {
            if (imuOdomQueue.front().header.stamp.toSec() < timeScanCur - 0.01)
                imuOdomQueue.pop_front();
            else
                break;
        }

        if (imuOdomQueue.empty())
        {
            ROS_INFO("empty Imu");
            return;
        }

        if (imuOdomQueue.front().header.stamp.toSec() > timeScanCur)
        {
            ROS_INFO("Late_Imu");
            return;
        }

        {
            // get start odometry at the beinning of the scan

            for (int i = 0; i < (int)imuOdomQueue.size(); ++i)
            {
                startImuOdomMsg = imuOdomQueue[i];

                if (ROS_TIME(&startImuOdomMsg) < timeScanCur)
                    continue;
                else
                    break;
            }
        }
        cloudInfo.imuAvailable = true;
        imuDeskewFlag = false;
        tf::quaternionMsgToTF(startImuOdomMsg.pose.pose.orientation, orientation);

        tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);

        odomImu[0] = startImuOdomMsg.pose.pose.position.x;
        odomImu[1] = startImuOdomMsg.pose.pose.position.y;
        odomImu[2] = startImuOdomMsg.pose.pose.position.z;
        odomImu[3] = roll;
        odomImu[4] = pitch;
        odomImu[5] = yaw;
        //? add: 同理，如果vins odom的增量平移变换不可用，则寻找imu odom的增量平移变换        
        //if (odomDeskewFlag == false)
        {
            // get end odometry at the end of the scan
            if (imuOdomQueue.back().header.stamp.toSec() < timeScanEnd)
                return;

            nav_msgs::Odometry endOdomMsg;

            for (int i = 0; i < (int)imuOdomQueue.size(); ++i)
            {
                endOdomMsg = imuOdomQueue[i];

                if (ROS_TIME(&endOdomMsg) < timeScanEnd)
                    continue;
                else
                    break;
            }

            //; 要保证前后的imu odom数据的id是一样的，即imu odom没有经过复位操作
            if (int(round(startImuOdomMsg.pose.covariance[0])) != int(round(endImuOdomMsg.pose.covariance[0])))
                return;

            Eigen::Affine3f transBegin = pcl::getTransformation(startImuOdomMsg.pose.pose.position.x, startImuOdomMsg.pose.pose.position.y, startImuOdomMsg.pose.pose.position.z, roll, pitch, yaw);
            tf::quaternionMsgToTF(endOdomMsg.pose.pose.orientation, orientation);
            
            tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);
            Eigen::Affine3f transEnd = pcl::getTransformation(endImuOdomMsg.pose.pose.position.x, endImuOdomMsg.pose.pose.position.y, endImuOdomMsg.pose.pose.position.z, roll, pitch, yaw);

            Eigen::Affine3f transBt = transBegin.inverse() * transEnd;

            pcl::getTranslationAndEulerAngles(transBt, odomIncreImu[0], odomIncreImu[1], odomIncreImu[2], odomIncreImu[3], odomIncreImu[4], odomIncreImu[5]);
        }
        odomDeskewFlag = true;
        imuDeskewFlag = true;
    }

    void findRotation(double pointTime, float *rotXCur, float *rotYCur, float *rotZCur)//AKF
    {
        *rotXCur = 0;
        *rotYCur = 0;
        *rotZCur = 0;

        int imuPointerFront = 0;
        while (imuPointerFront < imuPointerCur)
        {
            if (pointTime < imuTime[imuPointerFront])
                break;
            ++imuPointerFront;
        }

        if (pointTime > imuTime[imuPointerFront] || imuPointerFront == 0)
        {
            *rotXCur = imuRotX[imuPointerFront];
            *rotYCur = imuRotY[imuPointerFront];
            *rotZCur = imuRotZ[imuPointerFront];
        }
        else
        {
            int imuPointerBack = imuPointerFront - 1;
            double ratioFront = (pointTime - imuTime[imuPointerBack]) / (imuTime[imuPointerFront] - imuTime[imuPointerBack]);
            double ratioBack = (imuTime[imuPointerFront] - pointTime) / (imuTime[imuPointerFront] - imuTime[imuPointerBack]);
            *rotXCur = imuRotX[imuPointerFront] * ratioFront + imuRotX[imuPointerBack] * ratioBack;
            *rotYCur = imuRotY[imuPointerFront] * ratioFront + imuRotY[imuPointerBack] * ratioBack;
            *rotZCur = imuRotZ[imuPointerFront] * ratioFront + imuRotZ[imuPointerBack] * ratioBack;
        }
    }

    void findPosition(double relTime, float *posXCur, float *posYCur, float *posZCur)//AKF here
    {
        *posXCur = 0;
        *posYCur = 0;
        *posZCur = 0;

        // If the sensor moves relatively slow, like walking speed, positional deskew seems to have little benefits. Thus code below is commented.

        //? add: 打开去平移畸变，因为对于自动驾驶场景来说，高速状态下平移还是比较大的
        if(transDeskew)
        {
            if (cloudInfo.vinsOdomAvailable == false || cloudInfo.imuOdomAvailable == false || odomDeskewFlag == false)
                return;

            float ratio = relTime / (timeScanEnd - timeScanCur);

            *posXCur = ratio * odomIncreX;
            *posYCur = ratio * odomIncreY;
            *posZCur = ratio * odomIncreZ;
        }
    }

    PointType deskewPoint(PointType *point, double relTime)
    {
        if (deskewFlag == -1 || cloudInfo.imuAvailable == false)
            return *point;

        double pointTime = timeScanCur + relTime;

        float rotXCur, rotYCur, rotZCur;
        findRotation(pointTime, &rotXCur, &rotYCur, &rotZCur);

        float posXCur, posYCur, posZCur;
        findPosition(relTime, &posXCur, &posYCur, &posZCur);

        if (firstPointFlag == true)
        {
            transStartInverse = (pcl::getTransformation(posXCur, posYCur, posZCur, rotXCur, rotYCur, rotZCur)).inverse();
            firstPointFlag = false;
        }

        // transform points to start
        Eigen::Affine3f transFinal = pcl::getTransformation(posXCur, posYCur, posZCur, rotXCur, rotYCur, rotZCur);
        Eigen::Affine3f transBt = transStartInverse * transFinal;

        PointType newPoint;
        newPoint.x = transBt(0, 0) * point->x + transBt(0, 1) * point->y + transBt(0, 2) * point->z + transBt(0, 3);
        newPoint.y = transBt(1, 0) * point->x + transBt(1, 1) * point->y + transBt(1, 2) * point->z + transBt(1, 3);
        newPoint.z = transBt(2, 0) * point->x + transBt(2, 1) * point->y + transBt(2, 2) * point->z + transBt(2, 3);
        newPoint.intensity = point->intensity;

        return newPoint;
    }

    void projectPointCloud()
    {
        int cloudSize = laserCloudIn->points.size();
        // range image projection
        for (int i = 0; i < cloudSize; ++i)
        {
            PointType thisPoint;
            thisPoint.x = laserCloudIn->points[i].x;
            thisPoint.y = laserCloudIn->points[i].y;
            thisPoint.z = laserCloudIn->points[i].z;
            thisPoint.intensity = laserCloudIn->points[i].intensity;

            float range = pointDistance(thisPoint);
            if (range < lidarMinRange || range > lidarMaxRange)
                continue;

            int rowIdn = laserCloudIn->points[i].ring;
            if (rowIdn < 0 || rowIdn >= N_SCAN)
                continue;

            if (rowIdn % downsampleRate != 0)
                continue;

            int columnIdn = -1;
            if (sensor == SensorType::VELODYNE || sensor == SensorType::OUSTER)
            {
                float horizonAngle = atan2(thisPoint.x, thisPoint.y) * 180 / M_PI;
                static float ang_res_x = 360.0 / float(Horizon_SCAN);
                columnIdn = -round((horizonAngle - 90.0) / ang_res_x) + Horizon_SCAN / 2;
                if (columnIdn >= Horizon_SCAN)
                    columnIdn -= Horizon_SCAN;
            }
            else if (sensor == SensorType::LIVOX)
            {
                columnIdn = columnIdnCountVec[rowIdn];
                columnIdnCountVec[rowIdn] += 1;
            }

            if (columnIdn < 0 || columnIdn >= Horizon_SCAN)
                continue;

            if (rangeMat.at<float>(rowIdn, columnIdn) != FLT_MAX)
                continue;

            thisPoint = deskewPoint(&thisPoint, laserCloudIn->points[i].time);

            rangeMat.at<float>(rowIdn, columnIdn) = range;

            int index = columnIdn + rowIdn * Horizon_SCAN;
            fullCloud->points[index] = thisPoint;
        }
    }

    void cloudExtraction()
    {
        int count = 0;
        // extract segmented cloud for lidar odometry
        for (int i = 0; i < N_SCAN; ++i)
        {
            cloudInfo.startRingIndex[i] = count - 1 + 5;

            for (int j = 0; j < Horizon_SCAN; ++j)
            {
                if (rangeMat.at<float>(i, j) != FLT_MAX)
                {
                    // mark the points' column index for marking occlusion later
                    cloudInfo.pointColInd[count] = j;
                    // save range info
                    cloudInfo.pointRange[count] = rangeMat.at<float>(i, j);
                    // save extracted cloud
                    extractedCloud->push_back(fullCloud->points[j + i * Horizon_SCAN]);
                    // size of extracted cloud
                    ++count;
                }
            }
            cloudInfo.endRingIndex[i] = count - 1 - 5;
        }
    }

    void publishClouds()
    {
        cloudInfo.header = cloudHeader;
        cloudInfo.cloud_deskewed = publishCloud(pubExtractedCloud, extractedCloud, cloudHeader.stamp, lidarFrame);
        pubLaserCloudInfo.publish(cloudInfo);
    }
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "lvi_sam");

    ImageProjection IP;

    ROS_INFO("\033[1;32m----> Lidar Cloud Deskew Started.\033[0m");

    ros::MultiThreadedSpinner spinner(3);
    spinner.spin();

    return 0;
}