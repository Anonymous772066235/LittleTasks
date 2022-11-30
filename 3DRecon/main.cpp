#include"ctime"
#include "PointCloud2DEM.h"

int main()
{
        clock_t startTime;
        clock_t endTime;

        startTime=clock();
        std::cout<<"begin\n";
        //std::string PointCloudPath = "groundPointCloud_test.txt";
        std::string GroundPointCloudPath = "/home/wj/STUDY/ngrest_project/Proj_3D_XG/groundPointCloud_test.txt";
        std::string OutputDEMPath = "/home/wj/STUDY/ngrest_project/Proj_3D_XG/DEM_retest.tif";
        
        //PointCloud2DEM::Points2GroundPoints(PointCloudPath, GroundPointCloudPath, true, 0.5, 0.5);
        std::cout<<"16\n";
        PointCloud2DEM::GroundPoints2DEM(GroundPointCloudPath, OutputDEMPath,1);

        endTime=clock();

        std::cout<<(double)(endTime-startTime)/CLOCKS_PER_SEC<<std::endl;
        

        return 0;
}
