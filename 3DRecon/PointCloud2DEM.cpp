#include "PointCloud2DEM.h"
#include "BasicFunction.h"
#include "COORCONV.h"

void PointCloud2DEM::Points2GroundPoints(std::string PointCloudPath, std::string projectPath, bool CSF_SloopSmooth, double CSF_ClothResolution, double CSF_ClassThreshold)
{
        CSF csf;
        //读取txt点云文件
        csf.setProjPath(projectPath);
        csf.readPointsFromFile(PointCloudPath);
        // csf参数设置
        csf.params.bSloopSmooth = CSF_SloopSmooth;
        csf.params.cloth_resolution = CSF_ClothResolution;
        csf.params.class_threshold = CSF_ClassThreshold;
        csf.params.rigidness = 3;
        csf.params.interations = 500;
        csf.params.time_step = 0.65;
        //设置地面点与非地面点索引
        std::vector<int> GroundPointsIndex, OffGroundPointsIndex;
        //执行CSF
        csf.do_filtering(GroundPointsIndex, OffGroundPointsIndex, false);
        //保存地面点
        csf.savePoints(GroundPointsIndex, projectPath+"/GroundPointCloud.txt");
}

void PointCloud2DEM::GroundPoints2DEM(std::string GroundPointCloudPath,std::string projectPath,  double GridWidth)
{

        const int MAX_LENGTH = 512;
        char buffer[MAX_LENGTH];
        getcwd(buffer, 512);
        std::string current_path = buffer;

        BasicFunction::Write2Text(current_path+"/DEM_Path.txt",projectPath+"/Dem.tif");

        //从地面点云txt文件中按行读取xyz坐标至 std::vector<cv::Point3f> GroundPointCloud
        std::vector<cv::Point3f> GroundPointCloud; //建立一个点云存储器
        std::string aPointLine;
        cv::Point3f temp;
        std::ifstream PointsFile;

        PointsFile.open(GroundPointCloudPath);  //打开文件
        while (getline(PointsFile, aPointLine)) //按行读取txt文件
        {
                std::stringstream aLIneXYZ(aPointLine);
                // aLIneXYZ<<std::fixed<<std::setprecision(12)<<aPointLine;
                aLIneXYZ >> temp.x;
                aLIneXYZ >> temp.y;
                aLIneXYZ >> temp.z;

                GroundPointCloud.push_back(temp);
                aLIneXYZ.str("");
        }
        PointsFile.close();

        //日志进度
        std::string RunningStateMsg=BasicFunction::Convert2Json_new(200,"DEM_Interpolation is running",0.1,"false");
        BasicFunction::Write2Text(projectPath+"/RunningState_Log.txt",RunningStateMsg);


        //换3通道为单通道
        cv::Mat pointsMatXYZ = cv::Mat(GroundPointCloud).reshape(1);
        pointsMatXYZ.convertTo(pointsMatXYZ, CV_32F);

        //找点云坐标XY的最小值和最大值
        cv::Mat pointsMatX = pointsMatXYZ.col(0);
        cv::Mat pointsMatY = pointsMatXYZ.col(1);
        double x_min, y_min, x_max, y_max;
        cv::minMaxLoc(pointsMatX, &x_min, &x_max);
        cv::minMaxLoc(pointsMatY, &y_min, &y_max);

        //以XY最小值为原点，求点云相对原点的XY坐标
        for (int i = 0; i < pointsMatXYZ.rows; i++)
        {
                float *MatPtr = pointsMatXYZ.ptr<float>(i);
                MatPtr[0] -= x_min;
                MatPtr[1] -= y_min;
        }
        cv::Mat PointsGridXY = pointsMatXYZ(cv::Range::all(), cv::Range(0, 2));

        //日志进度
        RunningStateMsg=BasicFunction::Convert2Json_new(200,"DEM_Interpolation is running",0.3,"false");
        BasicFunction::Write2Text(projectPath+"/RunningState_Log.txt",RunningStateMsg);

        //计算二维格网行列数
        int ROWS = (y_max - y_min) / GridWidth + 1;
        int COLS = (x_max - x_min) / GridWidth + 1;

        // DEM存储器
        cv::Mat MatDEM(ROWS, COLS, CV_32F);

        //建立格网化点云的kdtree
        cv::flann::KDTreeIndexParams IndexParams(2);
        cv::flann::Index kdtree(PointsGridXY, IndexParams);
        unsigned queryNum = 3;
        std::vector<float> vecQuery(2);
        std::vector<int> vecIndex(queryNum);
        std::vector<float> vecDist(queryNum);
        cv::flann::SearchParams params(32);

        //日志进度
        RunningStateMsg=BasicFunction::Convert2Json_new(200,"DEM_Interpolation is running",0.5,"false");
        BasicFunction::Write2Text(projectPath+"/RunningState_Log.txt",RunningStateMsg);

        //基于反距离加权实现DEM插值
        std::vector<float> GridZ;
        float tempZ;
        float dy = GridWidth / 2;

        for (int i = 0; i < ROWS; i++)
        {
                float dx = GridWidth / 2;
                for (int j = 0; j < COLS; j++)
                {
                        vecQuery[0] = (float)dx;
                        vecQuery[1] = (float)dy;
                        kdtree.knnSearch(vecQuery, vecIndex, vecDist, queryNum, params);
                        tempZ = (pointsMatXYZ.ptr<float>(vecIndex[0])[2] / (vecDist[0]) + pointsMatXYZ.ptr<float>(vecIndex[1])[2] / (vecDist[1]) + pointsMatXYZ.ptr<float>(vecIndex[2])[2] / (vecDist[2])) / (1 / (vecDist[0]) + 1 / (vecDist[1]) + 1 / (vecDist[2]));
                        // tempZ = (pointsMatXYZ.at<float>(vecIndex[0], 2) / (vecDist[0]) + pointsMatXYZ.at<float>(vecIndex[1], 2) / (vecDist[1]) + pointsMatXYZ.at<float>(vecIndex[2], 2) / (vecDist[2])) / (1 / (vecDist[0]) + 1 / (vecDist[1]) + 1 / (vecDist[2]));
                        GridZ.push_back(tempZ);
                        dx += GridWidth;
                }
                dy += GridWidth;
        }


        //日志进度
        RunningStateMsg=BasicFunction::Convert2Json_new(200,"DEM_Interpolation is running",0.7,"false");
        BasicFunction::Write2Text(projectPath+"/RunningState_Log.txt",RunningStateMsg);

        //赋值给DEM存储器
        int k = 0;
        for (int i = ROWS-1; i >0; i--)
        {
                float *MatDEMPtr = MatDEM.ptr<float>(i);
                for (int j = 0; j < COLS; j++)
                {
                        MatDEMPtr[j] = GridZ[k];
                        k++;
                }
        }
        MatDEM.convertTo(MatDEM, CV_32F);
        //输出MatDEM为TIFF文件
        double GeoTrans[6] = {0, GridWidth, 0, 0, 0, GridWidth};

        
        std::ifstream fin(current_path+"/PointCloud_Path.txt");
        fin >> buffer;
        fin.close(); 


        std::vector<std::string>  Path_tmp;
        BasicFunction::Stringsplit(buffer,"/PointCloud.txt",Path_tmp);
        std::string pre_Path=Path_tmp[0];

        std::ifstream fin2(pre_Path+"/CornerPoints.txt");
        fin2 >> buffer;
        fin2.close(); 

        int average_lon=0;
        // 3.以"\""分割整个字符串为vector

        std::vector<std::string>  CornerVec;
        BasicFunction::Stringsplit(buffer,"\"",CornerVec);
        average_lon+=std::stof(CornerVec[1]);
        average_lon+=std::stof(CornerVec[5]);
        average_lon+=std::stof(CornerVec[9]);
        average_lon+=std::stof(CornerVec[13]);
        average_lon=average_lon/4.0;

        int zone=floor(average_lon/6)+31;

        // 把左上角的点转回经纬度
        WGS84Corr lonlat_TopLeft;
	UTMXYToLatLon(x_min, y_max, zone,false,lonlat_TopLeft);
	double lat_TopLeft = RadToDeg(lonlat_TopLeft.lat);  //注意弧度转角度！！！
	double lon_TopLeft = RadToDeg(lonlat_TopLeft.log);

        
        std::cout<<"\nlon_TopLeft:\t "<< lon_TopLeft<<std::endl;
        std::cout<<"lat_TopLeft:\t "<< lat_TopLeft<<std::endl;

        // 把右下角的点转回经纬度
        WGS84Corr lonlat_BottomRight;
	UTMXYToLatLon(x_max, y_min, zone,false,lonlat_BottomRight);
	double lat_BottomRight = RadToDeg(lonlat_BottomRight.lat);  //注意弧度转角度！！！
	double lon_BottomRight = RadToDeg(lonlat_BottomRight.log);

        std::cout<<"\nlon_BottomRight:\t "<< lon_BottomRight<<std::endl;
        std::cout<<"lat_BottomRight:\t "<< lat_BottomRight<<std::endl;


        // 写到6参数里面
        GeoTrans[0]=lon_TopLeft;
        GeoTrans[1]=(lon_BottomRight-lon_TopLeft)/ double(COLS);
        GeoTrans[3]=lat_TopLeft;
        GeoTrans[5]=(lat_BottomRight-lat_TopLeft)/ double(ROWS);




        OGRSpatialReference *oSRS = new OGRSpatialReference();
        oSRS->SetWellKnownGeogCS("WGS84");
        char *Proj = new char[1000];
        oSRS->exportToWkt(&Proj);

        std::vector<cv::Mat> VMatDEM;
        VMatDEM.push_back(MatDEM);
        ImageAndGeoinfo DEM;
        DEM._Image = VMatDEM;
        DEM._GeoTrans = GeoTrans;
        DEM._Proj = Proj;
        BasicFunction::SaveImageAsTiff(DEM._Image, projectPath+"/Dem.tif", DEM._Proj, DEM._GeoTrans);



        //日志进度
        RunningStateMsg=BasicFunction::Convert2Json_new(200,"DEM_Interpolation is running",0.9,"false");
        BasicFunction::Write2Text(projectPath+"/RunningState_Log.txt",RunningStateMsg);

        std::cout<<"The DEM file is:\n";
        std::cout<<projectPath+"/Dem.tif"<<std::endl;



        time_t timep;
        struct tm* p;
        time(&timep); //获取从1970至今过了多少秒，存入time_t类型的timep
        p = localtime(&timep);//用localtime将秒数转化为struct tm结构体

        //将时间转string
        std::string sTime = std::to_string(p->tm_year + 1900) + "-" + std::to_string(p->tm_mon+1) + "-" + std::to_string(p->tm_mday) +
                " " + std::to_string(p->tm_hour) + ":" + std::to_string(p->tm_min) + ":" + std::to_string(p->tm_sec);

        std::vector<std::string> pathsplit;
        BasicFunction::Stringsplit(projectPath+"/Dem.tif",".",pathsplit);
        std::string txtpath="";
        for(int k=0;k<pathsplit.size()-1;k++)
        {
              txtpath=txtpath+  pathsplit[k];
              if(k>0)
              {
                txtpath=txtpath+ ".";
              }
        }
        std::cout<<"The Auxiliary file is:\n";
        std::cout<<txtpath+".json"<<std::endl;
        std::ofstream outFile;             // output the Auxiliary file
        outFile.open(txtpath+".json");  // associate with a file
        outFile<<"{\"Data\":\"SZGCMX\",\"Format\":\".tif\",\"Name\":\"Dem.tif\",\"Time\":,\""<<sTime<<"\",\"Resolution\":\""<<GeoTrans[1]<<"*"<<GeoTrans[5]<<"\"}";
        outFile.close();

}
