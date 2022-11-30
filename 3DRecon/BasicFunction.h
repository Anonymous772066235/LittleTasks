#pragma once


#include <iostream>
#include <stdlib.h> 
#include <io.h>
#include <vector>
#include <fstream>
//#ifdef WIN32
// #include <direct.h>
//#endif
//UBUNTU 下的io库



#include <gdal.h>
#include <gdal_priv.h>
#include <gdal_alg.h>
#include <gdal_mdreader.h>
#include <cpl_conv.h>
#include <cpl_string.h>

#include "cpl_multiproc.h"
//#include "commonutils.h"
//#include "gdal_utils_priv.h"


#include <opencv2/core/core.hpp>
#include "opencv2/highgui/highgui.hpp"
#include <opencv2/features2d.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/cudafeatures2d.hpp>
#include <opencv2/xfeatures2d/cuda.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/xfeatures2d.hpp>
#include <opencv2/calib3d/calib3d.hpp>
//#include <opencv2/calib3d.hpp>
//#include "gms_matcher.h"
#include "tinyxml2.h"

//  for  readXML4Boundary
//  sudo apt-get install libxml2 
//  sudo ln -s /usr/include/libxml2/libxml   /usr/include/libxml
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/tree.h>




struct SATPoint2D
{
	// ------>sample
	// |
	// |
	// |
	// line

	double sample, line;
};

struct SATPoint3D
{
	//
	double L, P, H;
};

typedef struct RPCcoeffcient
{
	double LINE_OFF;
	double SAMP_OFF;
	double LAT_OFF;
	double LONG_OFF;
	double HEIGHT_OFF;
	double LINE_SCALE;
	double SAMP_SCALE;
	double LAT_SCALE;
	double LONG_SCALE;
	double HEIGHT_SCALE;
	double LINE_NUM_COEFF[20];
	double LINE_DEN_COEFF[20];
	double SAMP_NUM_COEFF[20];
	double SAMP_DEN_COEFF[20];
};


//利用RPC 投影，还要使用以下仿射变换变换到影像空间才能获取对应像点的坐标
struct RPCImAffine
{
	//L=lineb0+lineb1*L+lineb2*S
	//S=samplea0+samplea1*S+samplea2*L
	double samplea0; //shift parameters
	double samplea1;
	double samplea2;

	double lineb0; //shift parameters
	double lineb1;
	double lineb2;
};


//影像点结构体
struct ImagePoint
{
	float X, Y;//影像坐标
	short ID;//影像编号
	ImagePoint(float x, float y, short id)
	{
		X = x;
		Y = y;
		ID = id;
	}
};

struct FloatSort
{
	int index;
	float value;
};

struct ImageAndGeoinfo
{
	std::vector<cv::Mat> _Image;
	char * _Proj;
	double * _GeoTrans;
};

//xml中重要的元数据
struct GF2Metadata
{
	std::string SatelliteID;
	std::string SensorID;
	std::string ReceiveTime;
	std::string ProductFormat;
	int Bands;
	float ImageGSD;
	float RollViewingAngle;
	float PitchViewingAngle;
	float RollSatelliteAngle;
	float PitchSatelliteAngle;
	float YawSatelliteAngle;
	float SolarAzimuth;
	float SolarZenith;
	float SatelliteAzimuth;
	float SatelliteZenith;
	float IntegrationTime;
	std::string IntegrationLevel;
	std::string EarthEllipsoid;
	std::string MtfCorrection;
	float TopLeftLatitude;
	float TopLeftLongitude;
	float TopRightLatitude;
	float TopRightLongitude;
	float BottomRightLatitude;
	float BottomRightLongitude;
	float BottomLeftLatitude;
	float BottomLeftLongitude;

};
//影像边界
struct Boundary
{
	float TopLeftLatitude;
	float TopLeftLongitude;
	float TopRightLatitude;
	float TopRightLongitude;
	float BottomRightLatitude;
	float BottomRightLongitude;
	float BottomLeftLatitude;
	float BottomLeftLongitude;
};
//dem边界
struct DEMBoundary
{
	double MaxP, MaxL, MaxH;
	double MinP, MinL, MinH;
};


namespace BasicFunction
{

	//--------------------------------------------------------------------------------------------功能----------------------------------------------------------------------------------------------
	//GDAL格式转cv::Mat
	bool GDAL2Mat(GDALDataset* InputImageData, std::vector<cv::Mat>& OutputImage, char* Proj, double* GeoTransform);
	//Mat2GDAL
	GDALDataset* Mat2GDAL(std::vector<cv::Mat>& IntputImage, char* Proj, double* GeoTransform);

	//从影像文件中获取xml文件名
	std::string GetXmlFilePath(std::string ImageFilePath);
	//读取xml文件
	bool LoadXmlFile(std::string XmlFilePath, GF2Metadata& MetaData);

	//从metadata读取边界
	Boundary GetBoundaryFromMetadata(GF2Metadata Metadata);
	//获取影像边界
	DEMBoundary GetDEMBoundary(Boundary pImageBoundary);
	//从GeotransformParams读取边界
	Boundary GetBoundaryFromGeotransformParams(double * GeotransformParams, cv::Size ImageSize);
	//从多项式变换中获取边界
	DEMBoundary GetBoundaryFromPolynomialParams(cv::Size ImageSize, cv::Mat PolynomialParams);

	//constchar转char
	void ConstChar2Char(const char* Input, char* Output);

	//多项式变换
	cv::Point2f PolynomialTransform(float x, float y, cv::Mat PolynomialParams);

	//计算二维欧式距离
	double CalEuDistance2D(double x1, double y1, double x2, double y2);
	//特征提取
	void RobotSURF_CUDAMatch(const cv::Mat img1, const cv::Mat img2, std::vector<cv::Point2f>& p1, std::vector<cv::Point2f>& p2);
	//特征提取
	void RobotSURF_CUDAMatch(const cv::Mat img1, const cv::Mat img2, std::vector<cv::Point2f>& p1,
		std::vector<cv::Point2f>& p2, int patchheight, int patchwidth);
	//影像直方图匹配-------------------------------------------------------------
	//cv::Mat CHistMatch(cv::Mat img1, cv::Mat img2, int code = 8);
	cv::Mat CHistMatch(const cv::Mat img1, const cv::Mat img2);
	bool compare(FloatSort a, FloatSort b);
	int GetMinNumIndex(std::vector<float> data);
	//std::vector<float> CInterP(std::vector<float> data1, std::vector<float> data2);

	//分开文件夹路径和后缀名
	bool SegFilePath(std::string FilePath, std::string &FolderPath, std::string &FileName);
	//读取指定后缀名的文件
	void LoadFilesPath(std::string FolderPath, std::string Exd, std::vector<std::string>& FilesPath);
	//更改文件后缀名
	bool ChangeSuffix(std::string FilePath, std::string Suffix, std::string& NFilePath);
	//读取文件后缀名
	std::string GetFileSuffix(std::string FilePath);
	//按行读取txt中路径
	bool ReadImageNames(std::string ListPath, std::vector<std::string>& Data);
	//读取上级文件夹目录 
	bool GetFolderPath(std::string FilePath, std::string &FolderPath);
	//去除文件中一些空格
	void Trim(std::string & str);

	//--------------------------------------------------------------------------------------------读-----------------------------------------------------------------------------------------------
	//读取GTiff影像
	GDALDataset* ReadGTiffImage(std::string ImagePath, GDALRPCInfo& RPCParams);
	//读取多波段影像（多个文件）
	ImageAndGeoinfo ReadMulBandImageFromFolder(std::string ImagePath, GDALRPCInfo& RPCParams);
	//读取HDF格式文件
	GDALDataset* ReadHDFImage(std::string ImagePath);

	//bageyalu--------------------------------------------------------------------------------------
	//从xml读取边界
	bool GetBoundaryfromXML(std::string xmlpath, Boundary& bound);
	//读取单波段tif到mat
	cv::Mat ReadTifoneBand(std::string filename);
	//读取多波段tif到mat
	std::vector<cv::Mat> ReadTif(std::string filename);
	//投影系与坐标系转换-----------------------------------------------------------------
	cv::Point3f WGS842Proj(int des_code, double lon, double lat,double H);
	cv::Point2d Proj2WGS84(int code, cv::Point2d p);
	cv::Point2d WGS842Proj(std::string mode, double lon, double lat);//x、y
	cv::Point2d Proj2WGS84(std::string mode, cv::Point2d p, double mean_lon);//lon、lat
	int search_epsg_code(std::string datum,std::string mode, double lon);

	//字符串分割
	void Stringsplit(const std::string& str, const std::string& splits, std::vector<std::string>& res);

	//--------------------------------------------------------------------------------------------写-----------------------------------------------------------------------------------------------
	//将cv::Mat写入为tiff
	void SaveImageAsTiff(std::vector<cv::Mat> Image, std::string SaveFilePath, char* Proj, double* GeoTransform);
	//以指定格式写出影像，除开JPEG2000
	void SaveImage(std::vector<cv::Mat> Image, std::string Format, std::string SaveFilePath, char* Proj, double* GeoTransform);
	//以JPEG2000写出
	void SaveImageASOtherFormat(std::vector<cv::Mat> Image, std::string Format, std::string SaveFilePath, char* Proj, double* GeoTransform);
	//----------------------------------------------------------------------------------------------日志--------------------------------------------------------------------------------------------
	//创建文件夹
	int CreateFolder(std::string FolderPath);
	//删除文件夹
	int DeleteFolder(std::string FolderPath);
	//写入指定string到文件(没有则创建，不覆盖)
	void Write2Text(std::string FilePath, std::string Data);
	//删除指定文件
	int DeleteFile(std::string FilePath);
	//查询文件是否存在
	bool File_Exist(const std::string& name);
	//查询文件夹是否存在
	bool Path_Exist(const std::string& name);
	//将指定命令转化为json格式，并存到指定文件
	std::string Convert2Json_new(int Code, std::string Message, float Process, std::string isEnd);
	//original version of Convert2Json
	std::string Convert2Json(int Code, std::string Message, float Process = 0, std::string CurrentWork = "");
	//线性拉伸
	cv::Mat LinearStrech(cv::Mat img, float ratio=0.02, int code=8, int dtype= CV_8UC1);
	//读取底图块,前三个是rgb，后两个是x、y
	std::vector<cv::Mat> GetSubPatch(std::string Img_path, cv::Point2d TLPoint, cv::Point2d BRPoint);
	std::vector<cv::Mat> GetSubPatch(std::string Img_path, cv::Point2d TLPoint, int width, int height);
	//获取底图对应的坐标系统
	std::string GetBaseMapMode(std::string imgpath);
	//读取ZY3的xml中的角点
	int readXML4Boundary(std::string xmlPath, std::vector<std::vector<double>> & CornerPoints );
};

