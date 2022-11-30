#pragma once
//#include "stdafx.h"
//#include "RPCProcessing.h"
#include <vector>
#include <stdio.h>
#include <fstream>
#include <iostream>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc_c.h>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "BasicFunction.h"


using namespace std;



// typedef struct RPCCOEFFCIENT
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
// }RPCcoeffcient;


// // 利用RPC 投影，还要使用以下仿射变换变换到影像空间才能获取对应像点的坐标
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


// // 影像点结构体
// struct ImagePoint
// {
// 	float X, Y;//影像坐标
// 	short ID;//影像编号
// 	ImagePoint(float x, float y, short id)
// 	{
// 		X = x;
// 		Y = y;
// 		ID = id;
// 	}
// };


//将两张影像利用RPC参数进行生成核线影像，获取核线影像上同名点之后，将同名点坐标转换回去

class GetEpipolarImage
{
	
public:
	GetEpipolarImage();
	~GetEpipolarImage();
	void calCounterpartLines(SATPoint2D point, int layer);
	
	vector<vector<float>> createEpipolarImage(cv::Mat &le, cv::Mat &re,int type = 1);
	
	void calOriginalCoord();
	
	void RPCInterSection(SATPoint2D *lpt, SATPoint2D *rpt, int npt, SATPoint3D *ptObj);
	
	void setLRImage(cv::Mat limg, cv::Mat rimg);
	
	void setRPCcoefAndRPCImAffine(RPCcoeffcient lrpccoef, RPCcoeffcient rrpccoef);
	
	void createEpipolarImage1();
	
	void calOriginalCoord1(cv::Mat disp, vector<vector<float>>kb, cv::Mat &x, cv::Mat &y);
	
	void setoffset(double x, double y);

	RPCcoeffcient readrpc(char* filename, cv::Mat rpcdata);
private:
	float * leftImage;
	float * rightImage;
	double offx, offy;

	IplImage *leftImg, *rightImg;

	RPCcoeffcient lCoefficient;
	RPCcoeffcient rCoefficient;

	RPCImAffine lAffine;
	RPCImAffine rAffine;

	int lineCount = 0;
	int sampleCount = 0;

	double lk = 0, lb = 0;
	double rk = 0, rb = 0;
	cv::Mat lrgb, rrgb;

	vector<vector<float>> kb;
	unsigned char* rightEpipolarImage;
	unsigned char* leftEpipolarImage;
	void calCoefficients(vector<SATPoint2D> points, double& k, double& b);

	/*反算原图像坐标*/
	double ** lineKB;
	cv::Mat disp8;
	
	int inv(double *m1, int n);
	void matrixmulti(double *r1, double *r2, double *r, int m, int n, int p);
	double getaccumulation1(double rpc[20], SATPoint3D &ObjPt);
	double getaccumulation1(double *rpc, double L, double P, double H);
	bool GetInverseAffPara(RPCImAffine &srcAffPara, RPCImAffine &dstAffPara);
	double getpartialderivativeofL1(double *rpc, double L, double P, double H);
	double getpartialderivativeofP1(double *rpc, double L, double P, double H);
	double getpartialderivativeofH1(double *rpc, double L, double P, double H);

	void RPCImg2Obj(RPCcoeffcient& RPCcoef, double H, RPCImAffine& affpara, SATPoint2D pimgpt, SATPoint3D& ObjPt);
	SATPoint2D RPCObj2Img(RPCcoeffcient& RPCcoef, SATPoint3D& ObjPt, RPCImAffine& affpara);
	double GetPartialDerivativeofL(double Numrpc[20], double Denrpc[20], SATPoint3D& objpt, double SL);
	double GetPartialDerivativeofP(double Numrpc[20], double Denrpc[20], SATPoint3D& objpt, double SL);

	int ptype;


};

