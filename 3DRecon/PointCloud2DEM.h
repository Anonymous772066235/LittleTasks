#include "CSF.h"
#include "iostream"
#include "vector"
#include "iomanip"
#include "sstream"
#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/flann/miniflann.hpp"
#include "algorithm"
#include "cmath"
#include "BasicFunction.h"

namespace PointCloud2DEM
{

        void Points2GroundPoints(std::string PointCloudPath, std::string projectPath, bool CSF_SloopSmooth, double CSF_ClothResolution, double CSF_ClassThreshold);

        void GroundPoints2DEM(std::string GroundPointsPath,std::string projectPath,  double GridWidth);

}