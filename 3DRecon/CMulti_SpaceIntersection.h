#pragma once

#include <future>
#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <opencv2/xfeatures2d/nonfree.hpp>
#include "BasicFunction.h"

using namespace std;
// //
// //RPC参数
// // //RPCCOEFFCIENT
// struct RPCcoeffcient
// {
// 	double LINE_OFF;
// 	double SAMP_OFF;
// 	double LAT_OFF;
// 	double LONG_OFF;
// 	double HEIGHT_OFF;
// 	double LINE_SCALE;
// 	double SAMP_SCALE;
// 	double LAT_SCALE;
// 	double LONG_SCALE;
// 	double HEIGHT_SCALE;
// 	double LINE_NUM_COEFF[20];
// 	double LINE_DEN_COEFF[20];
// 	double SAMP_NUM_COEFF[20];
// 	double SAMP_DEN_COEFF[20];
// };

// //利用RPC投影，还要使用以下仿射变换变换到影像空间才能获取对应像点的坐标（优化）
// struct RPCImAffine
// {
// 	//L=lineb0+lineb1*L+lineb2*S
// 	//S=samplea0+samplea1*S+samplea2*L
// 	double samplea0; //shift parameters
// 	double samplea1;
// 	double samplea2;

// 	double lineb0; //shift parameters
// 	double lineb1;
// 	double lineb2;
// };

// //二维点结构体
// struct SATPoint2D
// {
// 	// ------>sample
// 	// |
// 	// |
// 	// |
// 	// line
// 	double sample, line;
// };
// //大地坐标
// struct SATPoint3D
// {
// 	//经纬度高程
// 	double L, P, H;
// };

// //影像点结构体
// struct ImagePoint
// {
// 	float X, Y;//影像坐标
// 	int ID;//影像编号
// 	ImagePoint(float x, float y, int id)
// 	{
// 		X = x;
// 		Y = y;
// 		ID = id;
// 	}
// };


class CMulti_SpaceIntersection
{

public:
	//构造、析构函数
	CMulti_SpaceIntersection();
	~CMulti_SpaceIntersection();
	//设置RPC参数
	void SetInputPathParams(std::vector<std::string>  ImageListPath);
	//设置其他参数
	void SetOtherParams(int T_ProjectErr, int k, double m,std::string projpath);
	//输入同名点
	void SetInputCorreImagePoints(std::vector<std::vector<ImagePoint>> CorrespondingImagePoints);
	//多片前方交会
	bool Multi_SpaceIntersection(std::vector<SATPoint3D> &PointCloud,std::string PointCloudSavePath="./PointClouds.txt");
private:
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	//路径参数

	// std::vector<std::string>  p_ImageListPath;
	std::vector<std::string> p_RPCParamsPaths;
	std::vector<std::string> p_RPCImAffineParamsPaths;


	//影像列表
	vector<std::string> p_ImageList;
	
	//左右影像RPC参数
	int p_T_ProjectErr;
	int p_k;
	double p_m;

	RPCcoeffcient p_lRPC, p_rRPC;
	RPCImAffine p_lAffine, p_rAffine;
	
	//总体RPC参数
	std::vector<RPCcoeffcient>p_RPCParams;
	std::vector<RPCImAffine>p_ImAffineParams;

	//工程目录
	std::string p_projpath;

	//同名点
	std::vector<std::vector<ImagePoint>> p_CorrespondingImagePoints;

	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	bool LoadRPCParams();
	//读
	std::vector<std::vector<double>> ReadRPCImAffineParamsFromTxt(std::string FilePath);
	RPCcoeffcient ReadRPC(char* filename);
	//分割
	void Stringsplit(const std::string& str, const std::string& splits, std::vector<std::string>& res);
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	//点云滤波
	vector<cv::Point3f> SOR(vector<cv::Point3f> InputPointCloud, int kdnum, double m);
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	//坐标转换
	SATPoint3D XYZ2BLH(double x, double y, double z);
	cv::Point3f BLH2XYZ(double b, double l, double h);
	void LonLat2UTM(double longitude, double latitude, double& UTME, double& UTMN);
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	//输入一个仿射变换参数 求出逆仿射变换参数
	bool GetInverseAffPara(RPCImAffine& srcAffPara, RPCImAffine& dstAffPara);

	//RPC参数
	void SetImgAffine(RPCImAffine lAffine, RPCImAffine rAffine);
	void SetRPCParams(RPCcoeffcient& lRPCcoef, RPCcoeffcient& rRPCcoef);

	//多片前方交会
	SATPoint3D Multi_RPCIntersection(vector<ImagePoint> sp,  int& flag);
	//前方交会
	bool RPCIntersection(SATPoint2D* lpt, SATPoint2D* rpt, int npt, SATPoint3D* ptObj);

	//一对影像计算A和L矩阵
	cv::Mat GetAandL(SATPoint2D* lpt, SATPoint3D* ptObj, int flag);
	//多张影像计算A和L
	cv::Mat Multi_GetAandL(vector<ImagePoint> sp, SATPoint3D* pt3d);
	//获取ARC
	cv::Mat GetAandRCandL(SATPoint2D* lpt, SATPoint3D* ptObj, RPCImAffine imaffine);

	//将地面点投影到影像
	SATPoint2D RPCObj2Img(RPCcoeffcient& RPCcoef, SATPoint3D& ObjPt, RPCImAffine& affpara);
	//给定地物点高程确定地面点坐标
	void RPCImg2Obj(RPCcoeffcient& RPCcoef, double H, RPCImAffine& affpara, SATPoint2D pimgpt, SATPoint3D& ObjPt);
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	//求逆
	int inv(double* m1, int n);
	//矩阵乘
	void matrixmulti(double* r1, double* r2, double* r, int m, int n, int p);

	//偏导数
	double getpartialderivativeofL1(double* rpc, double L, double P, double H);
	double getpartialderivativeofP1(double* rpc, double L, double P, double H);
	double getpartialderivativeofH1(double* rpc, double L, double P, double H);
	//求和
	double getaccumulation1(double* rpc, double L, double P, double H);
	double getaccumulation1(double rpc[20], SATPoint3D& ObjPt);
	//像方到物方
	double GetPartialDerivativeofP(double Numrpc[20], double Denrpc[20], SATPoint3D& objpt, double SL);
	double GetPartialDerivativeofL(double Numrpc[20], double Denrpc[20], SATPoint3D& objpt, double SL);


};

