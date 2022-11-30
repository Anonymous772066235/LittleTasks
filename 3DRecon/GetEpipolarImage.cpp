//#include "stdafx.h"
#include <string.h>
#include <stdio.h>
#include <cmath>
#include "GetEpipolarImage.h"
#include"CMulti_SpaceIntersection.h"


typedef unsigned char BYTE;
using namespace std;
using namespace cv;


GetEpipolarImage::GetEpipolarImage()
{
}
GetEpipolarImage::~GetEpipolarImage()
{
}


RPCcoeffcient GetEpipolarImage::readrpc(char* filename, cv::Mat rpcdata)
{
	RPCcoeffcient rpcdata1;
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
	long double* N_L;
	long double* D_L;
	long double* N_S;
	long double* D_S;
	FILE* fp = fopen(filename, "r");
	char tmp1[200], tmp2[200], tmp3[200];
	fscanf(fp, "%s %lf %s", tmp1, &LINE_OFF, tmp2);
	fscanf(fp, "%s %lf %s", tmp1, &SAMP_OFF, tmp2);
	fscanf(fp, "%s %lf %s", tmp1, &LAT_OFF, tmp2);
	fscanf(fp, "%s %lf %s", tmp1, &LONG_OFF, tmp2);
	fscanf(fp, "%s %lf %s", tmp1, &HEIGHT_OFF, tmp2);
	fscanf(fp, "%s %lf %s", tmp1, &LINE_SCALE, tmp2);
	fscanf(fp, "%s %lf %s", tmp1, &SAMP_SCALE, tmp2);
	fscanf(fp, "%s %lf %s", tmp1, &LAT_SCALE, tmp2);
	fscanf(fp, "%s %lf %s", tmp1, &LONG_SCALE, tmp2);
	fscanf(fp, "%s %lf %s", tmp1, &HEIGHT_SCALE, tmp2);
	N_L = new long double[20];
	D_L = new long double[20];
	N_S = new long double[20];
	D_S = new long double[20];
	for (int i = 0; i < 20; i++)
	{
		fscanf(fp, "%s %Lf", tmp1, &N_L[i]);

		//cout<< setprecision(44) <<N_L[i]<<endl;
	}
	for (int i = 0; i < 20; i++)
	{
		fscanf(fp, "%s %Lf", tmp1, &D_L[i]);
		//cout<<D_L[i]<<endl;
	}
	for (int i = 0; i < 20; i++)
	{
		fscanf(fp, "%s %Lf", tmp1, &N_S[i]);
		//cout<<N_S[i]<<endl;
	}
	for (int i = 0; i < 20; i++)
	{
		fscanf(fp, "%s %Lf", tmp1, &D_S[i]);
		//cout<<D_S[i]<<endl;
	}
	fclose(fp);

	Mat temp(90, 1, CV_32F);

	temp.at<float>(0, 0) = LINE_OFF;
	temp.at<float>(1, 0) = SAMP_OFF;
	temp.at<float>(2, 0) = LAT_OFF;
	temp.at<float>(3, 0) = LONG_OFF;
	temp.at<float>(4, 0) = HEIGHT_OFF;
	temp.at<float>(5, 0) = LINE_SCALE;
	temp.at<float>(6, 0) = SAMP_SCALE;
	temp.at<float>(7, 0) = LAT_SCALE;
	temp.at<float>(8, 0) = LONG_SCALE;
	temp.at<float>(9, 0) = HEIGHT_SCALE;

	for (int i = 0; i < 20; i++)
	{
		temp.at<float>(i + 10, 0) = N_L[i];
		temp.at<float>(i + 30, 0) = D_L[i];
		temp.at<float>(i + 50, 0) = N_S[i];
		temp.at<float>(i + 70, 0) = D_S[i];
		rpcdata1.LINE_DEN_COEFF[i] = D_L[i];
		rpcdata1.LINE_NUM_COEFF[i] = N_L[i];
		rpcdata1.SAMP_DEN_COEFF[i] = D_S[i];
		rpcdata1.SAMP_NUM_COEFF[i] = N_S[i];
	}
	rpcdata = temp;
	rpcdata1.HEIGHT_OFF = HEIGHT_OFF;
	rpcdata1.HEIGHT_SCALE = HEIGHT_SCALE;
	rpcdata1.LAT_OFF = LAT_OFF;
	rpcdata1.LAT_SCALE = LAT_SCALE;
	rpcdata1.LINE_OFF = LINE_OFF;
	rpcdata1.LINE_SCALE = LINE_SCALE;
	rpcdata1.LONG_OFF = LONG_OFF;
	rpcdata1.LONG_SCALE = LONG_SCALE;
	rpcdata1.SAMP_OFF = SAMP_OFF;
	rpcdata1.SAMP_SCALE = SAMP_SCALE;

	return rpcdata1;

}



//zys
//data 2012-12-12
// Fs = NumS(P,L,H)-Sample*DenS(P,L,H) = 0;
// Fl = NumL(P,L,H)-Line*DenL(P,L,H)=0;
//求与P相关的偏导数
double GetEpipolarImage::GetPartialDerivativeofP(double Numrpc[20], double Denrpc[20], SATPoint3D& objpt, double SL)
{
	double P = objpt.P;
	double L = objpt.L;
	double H = objpt.H;
	double sum = Numrpc[2] + Numrpc[4] * L + Numrpc[6] * H + 2 * Numrpc[8] * P +
		Numrpc[10] * L * H + 2 * Numrpc[12] * L * P + Numrpc[14] * L * L + 3 * Numrpc[15] * P * P + Numrpc[16] * H * H +
		2 * Numrpc[18] * P * H
		- SL * (Denrpc[2] + Denrpc[4] * L + Denrpc[6] * H + 2 * Denrpc[8] * P +
			Denrpc[10] * L * H + 2 * Denrpc[12] * L * P + Denrpc[14] * L * L + 3 * Denrpc[15] * P * P + Denrpc[16] * H * H +
			2 * Denrpc[18] * P * H);

	return sum;
}


//求与L相关的偏导数
double GetEpipolarImage::GetPartialDerivativeofL(double Numrpc[20], double Denrpc[20], SATPoint3D& objpt, double SL)
{
	double P = objpt.P;
	double L = objpt.L;
	double H = objpt.H;
	double sum = Numrpc[1] + Numrpc[4] * P + Numrpc[5] * H + 2 * Numrpc[7] * L +
		Numrpc[10] * P * H + 3 * Numrpc[11] * L * L + Numrpc[12] * P * P + Numrpc[13] * H * H + 2 * Numrpc[14] * L * P +
		2 * Numrpc[17] * L * H
		- SL * (Denrpc[1] + Denrpc[4] * P + Denrpc[5] * H + 2 * Denrpc[7] * L +
			Denrpc[10] * P * H + 3 * Denrpc[11] * L * L + Denrpc[12] * P * P + Denrpc[13] * H * H + 2 * Denrpc[14] * L * P +
			2 * Denrpc[17] * L * H);

	return sum;
}

//利用RPC参数将三维点投影到影像
SATPoint2D GetEpipolarImage::RPCObj2Img(RPCcoeffcient& RPCcoef, SATPoint3D& ObjPt, RPCImAffine& affpara)
{
	SATPoint2D ptImg;
	double P = (ObjPt.P - RPCcoef.LAT_OFF) / RPCcoef.LAT_SCALE;
	double L = (ObjPt.L - RPCcoef.LONG_OFF) / RPCcoef.LONG_SCALE;
	double H = (ObjPt.H - RPCcoef.HEIGHT_OFF) / RPCcoef.HEIGHT_SCALE;

	double M[20];
	M[0] = 1, M[1] = L, M[2] = P, M[3] = H,
		M[4] = L * P, M[5] = L * H, M[6] = P * H, M[7] = L * L, M[8] = P * P, M[9] = H * H,
		M[10] = P * L * H, M[11] = L * L * L, M[12] = L * P * P, M[13] = L * H * H, M[14] = L * L * P,
		M[15] = P * P * P, M[16] = P * H * H, M[17] = L * L * H, M[18] = P * P * H, M[19] = H * H * H;

	double NumL = 0, DenL = 0, NumS = 0, DenS = 0;
	for (int i = 0; i < 20; i++)
	{
		NumL += RPCcoef.LINE_NUM_COEFF[i] * M[i];
		DenL += RPCcoef.LINE_DEN_COEFF[i] * M[i];
		NumS += RPCcoef.SAMP_NUM_COEFF[i] * M[i];
		DenS += RPCcoef.SAMP_DEN_COEFF[i] * M[i];
	}

	double linetmp = (NumL / DenL) * RPCcoef.LINE_SCALE + RPCcoef.LINE_OFF;
	double sampetmp = (NumS / DenS) * RPCcoef.SAMP_SCALE + RPCcoef.SAMP_OFF;

	ptImg.sample = affpara.samplea0 + affpara.samplea1 * sampetmp + affpara.samplea2 * linetmp;
	ptImg.line = affpara.lineb0 + affpara.lineb1 * sampetmp + affpara.lineb2 * linetmp;
	return ptImg;//#
}


// 利用RPC将二维点转到三维（已知高程）
void GetEpipolarImage::RPCImg2Obj(RPCcoeffcient& RPCcoef, double H, RPCImAffine& affpara, SATPoint2D pimgpt, SATPoint3D& ObjPt)
{
	//坐标归一化

	//给定的仿射变换矩阵 的逆变换

	double sample = affpara.samplea0 + affpara.samplea1 * pimgpt.sample + affpara.samplea2 * pimgpt.line;
	double line = affpara.lineb0 + affpara.lineb1 * pimgpt.sample + affpara.lineb2 * pimgpt.line;

	sample = (sample - RPCcoef.SAMP_OFF) / RPCcoef.SAMP_SCALE;
	line = (line - RPCcoef.LINE_OFF) / RPCcoef.LINE_SCALE;

	ObjPt.L = 0;
	ObjPt.P = 0;

	ObjPt.H = (H - RPCcoef.HEIGHT_OFF) / RPCcoef.HEIGHT_SCALE;

	//

	//v=bx-lcoef
	double B[2][2], lcoef[2][1], x[2][1], BT[2][2], BTB[2][2], BTl[2][1];
	int noitrative = 0;

	//x[2][1], det P det L
	while (noitrative < 15)
	{
		B[0][0] = GetPartialDerivativeofL(RPCcoef.SAMP_NUM_COEFF, RPCcoef.SAMP_DEN_COEFF, ObjPt, sample);
		B[0][1] = GetPartialDerivativeofP(RPCcoef.SAMP_NUM_COEFF, RPCcoef.SAMP_DEN_COEFF, ObjPt, sample);
		B[1][0] = GetPartialDerivativeofL(RPCcoef.LINE_NUM_COEFF, RPCcoef.LINE_DEN_COEFF, ObjPt, line);
		B[1][1] = GetPartialDerivativeofP(RPCcoef.LINE_NUM_COEFF, RPCcoef.LINE_DEN_COEFF, ObjPt, line);

		lcoef[0][0] = -(getaccumulation1(RPCcoef.SAMP_NUM_COEFF, ObjPt) - sample * getaccumulation1(RPCcoef.SAMP_DEN_COEFF, ObjPt));
		lcoef[1][0] = -(getaccumulation1(RPCcoef.LINE_NUM_COEFF, ObjPt) - line * getaccumulation1(RPCcoef.LINE_DEN_COEFF, ObjPt));

		for (int p = 0; p < 2; p++)
		{
			for (int q = 0; q < 2; q++)
			{
				BT[q][p] = B[p][q];
			}
		}

		matrixmulti(*BT, *B, *BTB, 2, 2, 2);
		matrixmulti(*BT, *lcoef, *BTl, 2, 2, 1);


		if (inv(*BTB, 2))
		{
			matrixmulti(*BTB, *BTl, *x, 2, 2, 1);

		}
		else break;

		ObjPt.L += x[0][0];
		ObjPt.P += x[1][0];
		/*if(fabs(x[0][0])<1.0e-10 && fabs(x[1][0])<1.0e-10)
			break;*/
		noitrative++;
	}

	ObjPt.H = H;
	ObjPt.L = ObjPt.L * RPCcoef.LONG_SCALE + RPCcoef.LONG_OFF;
	ObjPt.P = ObjPt.P * RPCcoef.LAT_SCALE + RPCcoef.LAT_OFF;
}



//读取左右影像灰度值，并保存在float数组中
void GetEpipolarImage::setLRImage(Mat limg,Mat rimg)
{

	Mat lg, rg;
	if (limg.depth() == 3)
		cv::cvtColor(limg, lg, CV_BGR2GRAY);
	else
		lg = limg.clone();
	if (rimg.depth() == 3)
		cv::cvtColor(rimg, rg, CV_BGR2GRAY);
	else
		rg = rimg.clone();

	lrgb = limg;
	rrgb = rimg;

	int nImgSizeY1 = limg.rows;
	int nImgSizeX1 = limg.cols;

	leftImage = new float[nImgSizeX1*nImgSizeY1];//用一个float类型的数组存储图片的灰度值
	rightImage = new float[nImgSizeX1*nImgSizeY1];//用一个float类型的数组存储图片的灰度值
	//用一个float类型的数组存储图片的灰度值
	rightImage = new float[nImgSizeX1*nImgSizeY1];//用一个float类型的数组存储图片的灰度值
	lineCount = limg.rows;
	sampleCount = limg.cols;
	

	leftEpipolarImage = new unsigned char[nImgSizeX1*nImgSizeY1];
	rightEpipolarImage = new  unsigned char[nImgSizeX1*nImgSizeY1];


	for (int i = 0; i < sampleCount; i++) 
	{
		for (int j = 0; j < lineCount; j++) 
		{
			
			leftImage[sampleCount*j + i] = lg.at<uchar>(j,i);
			rightImage[sampleCount*j + i] = rg.at<uchar>(j, i);

			leftEpipolarImage[sampleCount*j + i] = 0;
			rightEpipolarImage[sampleCount*j + i] = 0;
		}
	}


}

//依据一系列points，拟合一个y=kx+b形式的直线
void GetEpipolarImage::calCoefficients(vector<SATPoint2D> points, double & k, double & b)
{
	double SigmaX(0);
	double SigmaY(0);
	double SigmaXY(0);
	double SigmaX2(0);
	int n = points.size();
	for (int i = 0; i < n; i++)
	{
		SigmaX += points[i].sample;
		SigmaY += points[i].line;
		SigmaXY += points[i].sample * points[i].line;
		SigmaX2 += points[i].sample * points[i].sample;
	}
	k = (n * SigmaXY - SigmaX * SigmaY) / (n * SigmaX2 - SigmaX * SigmaX);
	b = SigmaY / n - k * SigmaX / n;

	double k1 = (points[0].line - points[1].line) / (points[0].sample - points[1].sample);
	double b1 = points[0].line - k1*points[0].sample;

}

//读取相应的参数
void GetEpipolarImage::setRPCcoefAndRPCImAffine(RPCcoeffcient lrpccoef, RPCcoeffcient rrpccoef)
{
	lCoefficient=lrpccoef;
	rCoefficient=rrpccoef;

	lAffine.lineb0 = 0;
	lAffine.lineb2 = 1;
	lAffine.lineb1 = 0;
	lAffine.samplea0 = 0;
	lAffine.samplea2 = 0;
	lAffine.samplea1 = 1;

	rAffine.lineb0 = 0;
	rAffine.lineb2 = 1;
	rAffine.lineb1 = 0;
	rAffine.samplea0 = 0;
	rAffine.samplea2 = 0;
	rAffine.samplea1 = 1;

}


//计算某个point在右影像上的直线
//point:需要的点
//layer:网格的层数
void GetEpipolarImage::calCounterpartLines(SATPoint2D point, int layer)
{
	//RPCProcessing processing;
	vector<SATPoint2D> rPoints;

	double maxHeight = 2000;//设置高程阈值
	double minHeight = 0;
	for (int i = 0; i < layer; i++) 
	{
		double tmpHeight = i* (maxHeight - minHeight) / (layer*1.0) + minHeight;	//根据层数选取相应的高度H
		
		SATPoint3D tmp3Dpoint, tmp3Dpoint2;
		RPCImAffine inverseLAffine;	//三维算二维是此参数需要取逆
		GetInverseAffPara(lAffine, inverseLAffine);
		RPCImg2Obj(lCoefficient, tmpHeight, inverseLAffine, point, tmp3Dpoint);		//根据点的坐标和高度反算成实际的3维坐标
		SATPoint2D tmpRPoint = RPCObj2Img(rCoefficient, tmp3Dpoint, rAffine);	//三维点换算到右影像的二维点
		
		CvPoint tmpCvPoint;
		tmpCvPoint.x = tmpRPoint.sample;
		tmpCvPoint.y = tmpRPoint.line;

		rPoints.push_back(tmpRPoint);	//存入数组，以便计算直线
	}

	//计算右边对应点拟合的直线
	calCoefficients(rPoints, rk, rb);


	
	vector<SATPoint2D> lPoints;//左影像点集
	for (int i = 0; i < layer; i++)
	{
		SATPoint2D tmpRPoint2;//随机找一个右影像线上的点，这里是横坐标取全图的一半，纵坐标用y=kx+b计算
		tmpRPoint2.line = lineCount *i / (layer*1.0); //斜率太大，用y来求x
		tmpRPoint2.sample = tmpRPoint2.line / rk - rb / rk;

		SATPoint3D tmp3Dpoint;
		RPCImAffine inverseRAffine;
		GetInverseAffPara(rAffine, inverseRAffine);
		RPCImg2Obj(rCoefficient,200, inverseRAffine, tmpRPoint2, tmp3Dpoint);

		SATPoint2D tmpLPoint = RPCObj2Img(lCoefficient, tmp3Dpoint, lAffine); //左影像上的点
		lPoints.push_back(tmpLPoint);

		CvPoint tmpCvPoint;
		tmpCvPoint.x = tmpLPoint.sample;
		tmpCvPoint.y = tmpLPoint.line;

	}
	lPoints.push_back(point);

	//根据点结拟合成直线
	calCoefficients(lPoints, lk, lb);
}
//创建核线影像
vector<vector<float>> GetEpipolarImage::createEpipolarImage(Mat& le, Mat& re,int type )
{
	ptype = type;
	//用来显示核线影像

	char fileName[] = "../outputdata/kb/linekb.txt";
	ofstream fout(fileName);

	//创建核线影像容器
	le = Mat::zeros(lrgb.rows, lrgb.cols, CV_8UC1);
	re = Mat::zeros(rrgb.rows, rrgb.cols, CV_8UC1);
	
	double _min, _max;
	minMaxIdx(lrgb, &_min, &_max);

	//寻找核线影像对应点
	for (int j = 1; j < lrgb.cols; j++)
	{
		//样本点
		SATPoint2D itmp;
		itmp.sample = j;
		itmp.line = lrgb.rows / 2.0;

		//计算对应的核线
		calCounterpartLines(itmp, 10);

		//y=kx+b
		vector<float> kb_tmp;
		kb_tmp.push_back(lk);
		kb_tmp.push_back(lb);
		kb_tmp.push_back(rk);
		kb_tmp.push_back(rb);
		kb.push_back(kb_tmp);

		//写入到文件
		fout << j << " " << lk << " " << lb  << endl;

		//左影像纠正
		for (int i = 1; i <= lrgb.rows; i++)
		{
			//int lTmpY = lk * (tmpSample + offx) + lb - offy;
			int itmpx = i / lk - lb / lk;
			
			if (itmpx < lrgb.cols && itmpx > 0)
			{
				le.at<uchar>(i - 1, j - 1) = lrgb.at<uchar>(i - 1, itmpx - 1);
			}
			/*if (i == 2500 && j == 2500)
			{
				cout << i - 1 << "  " << itmpx - 1 << endl;
			}*/
		}

		//右影像纠正
		for (int i = 1; i <= rrgb.rows; i++)
		{
			//int lTmpY = lk * (tmpSample + offx) + lb - offy;
			int itmpx = i / rk - rb / rk;

			if ((itmpx < rrgb.cols && itmpx > 0))
			{
				re.at<uchar>(i - 1, j - 1) = rrgb.at<uchar>(i - 1, itmpx - 1);
			}
			/*if (i == 2500 && j == 2500)
			{
				cout << i - 1 <<"  " << itmpx - 1 << endl;
			}*/
		}

		
	}

	return kb;

}

//核线影像上的同名点坐标转换回去
void GetEpipolarImage::calOriginalCoord()
{
	//反算原图像坐标
	//测试时使用

	lineKB = new double*[sampleCount + 1];
	for (int i = 0; i < sampleCount + 1; i++)
		lineKB[i] = new double[4];
	

	SATPoint2D* lPoint2Ds = new SATPoint2D[sampleCount*lineCount];
	SATPoint2D* rPoint2Ds = new SATPoint2D[sampleCount*lineCount];
	for (int i = 0; i < sampleCount; i++) 
	{
		for (int j = 0; j < lineCount; j++) 
		{
			int lx = lineCount - 1 - i, ly = j;
			int rx = lineCount - 1 - i - disp8.at<float>(ly, lx), ry = ly;
			//for(int j=0;j<24;j++) {
			//	int lx = j*80, ly = 150;
			//	int rx = j*80+ disp8.at<float>(ly, lx), ry = ly;

			double centerX = sampleCount / 2.0;
			double centerY = lineCount / 2.0;
			int lx_90 = (ly - centerY) + centerX;
			int ly_90 = -(lx - centerX) + centerY;
			int rx_90 = (ry - centerY) + centerX;
			int ry_90 = -(rx - centerX) + centerY;

			CvPoint tmppoint;
			double lk = kb[lx_90][0];
			double lb = kb[lx_90][1];
			double rk = kb[lx_90][2];
			double rb = kb[lx_90][3];

			int originalLY = ly_90;
			int originalLX = ly_90 / lk - lb / lk;
			int originalRY = ry_90;
			int originalRX = ry_90 / rk - rb / rk;

			lPoint2Ds[i*lineCount + j].sample = originalLX;
			lPoint2Ds[i*lineCount + j].line = originalLY;

			rPoint2Ds[i*lineCount + j].sample = originalRX;
			rPoint2Ds[i*lineCount + j].line = originalRY;
		}
	}

	SATPoint3D* point3Ds = new SATPoint3D[sampleCount*lineCount];
	RPCInterSection(lPoint2Ds, rPoint2Ds, sampleCount*lineCount, point3Ds);
}

void GetEpipolarImage::RPCInterSection(SATPoint2D * lpt, SATPoint2D * rpt, int npt, SATPoint3D * ptObj)
{
	double A[4][3], l[4][1], d[3][1], AT[3][4], ATA[3][3], ATl[3][1];
	d[0][0] = 1;
	d[1][0] = 1;
	d[2][0] = 50;
	double Pn = 0, Ln = 0, Hn = 0;
	double NumL = 0, DenL = 0, NumS = 0, DenS = 0;
	double dNumLdLn = 0.0, dNumLdPn = 0.0, dNumLdHn = 0.0;
	double dDenLdLn = 0.0, dDenLdPn = 0.0, dDenLdHn = 0.0;
	double dNumSdLn = 0.0, dNumSdPn = 0.0, dNumSdHn = 0.0;
	double dDenSdLn = 0.0, dDenSdPn = 0.0, dDenSdHn = 0.0;  //上面四行用来记录正解法偏导数

	int numofiterative = 0;
	double Ltmp, Stmp;

	for (int i = 0; i<npt; i++)
	{
		ptObj[i].L = (lCoefficient.LONG_OFF + rCoefficient.LONG_OFF) / 2;
		ptObj[i].P = (lCoefficient.LAT_OFF + rCoefficient.LAT_OFF) / 2;
		ptObj[i].H = (lCoefficient.HEIGHT_OFF + rCoefficient.HEIGHT_OFF) / 2;
	}

	for (int i = 0; i<npt; i++)
	{
		// while(numofiterative<30&&(fabs(d[0][0])>0.00001||fabs(d[1][0])>0.00001||fabs(d[2][0])>0.001))
		while (numofiterative<15)
		{
			//……………………处理左片………………	

			Pn = (ptObj[i].P - lCoefficient.LAT_OFF) / lCoefficient.LAT_SCALE;
			Ln = (ptObj[i].L - lCoefficient.LONG_OFF) / lCoefficient.LONG_SCALE;
			Hn = (ptObj[i].H - lCoefficient.HEIGHT_OFF) / lCoefficient.HEIGHT_SCALE;

			NumL = getaccumulation1(lCoefficient.LINE_NUM_COEFF, Ln, Pn, Hn);
			DenL = getaccumulation1(lCoefficient.LINE_DEN_COEFF, Ln, Pn, Hn);
			NumS = getaccumulation1(lCoefficient.SAMP_NUM_COEFF, Ln, Pn, Hn);
			DenS = getaccumulation1(lCoefficient.SAMP_DEN_COEFF, Ln, Pn, Hn);

			dNumLdLn = getpartialderivativeofL1(lCoefficient.LINE_NUM_COEFF, Ln, Pn, Hn);
			dNumLdPn = getpartialderivativeofP1(lCoefficient.LINE_NUM_COEFF, Ln, Pn, Hn);
			dNumLdHn = getpartialderivativeofH1(lCoefficient.LINE_NUM_COEFF, Ln, Pn, Hn);

			dDenLdLn = getpartialderivativeofL1(lCoefficient.LINE_DEN_COEFF, Ln, Pn, Hn);
			dDenLdPn = getpartialderivativeofP1(lCoefficient.LINE_DEN_COEFF, Ln, Pn, Hn);
			dDenLdHn = getpartialderivativeofH1(lCoefficient.LINE_DEN_COEFF, Ln, Pn, Hn);

			dNumSdLn = getpartialderivativeofL1(lCoefficient.SAMP_NUM_COEFF, Ln, Pn, Hn);
			dNumSdPn = getpartialderivativeofP1(lCoefficient.SAMP_NUM_COEFF, Ln, Pn, Hn);
			dNumSdHn = getpartialderivativeofH1(lCoefficient.SAMP_NUM_COEFF, Ln, Pn, Hn);

			dDenSdLn = getpartialderivativeofL1(lCoefficient.SAMP_DEN_COEFF, Ln, Pn, Hn);
			dDenSdPn = getpartialderivativeofP1(lCoefficient.SAMP_DEN_COEFF, Ln, Pn, Hn);
			dDenSdHn = getpartialderivativeofH1(lCoefficient.SAMP_DEN_COEFF, Ln, Pn, Hn);

			A[0][0] = (lCoefficient.LINE_SCALE / lCoefficient.LONG_SCALE)*((dNumLdLn*DenL - dDenLdLn*NumL) / (DenL*DenL));
			A[0][1] = (lCoefficient.LINE_SCALE / lCoefficient.LAT_SCALE)*((dNumLdPn*DenL - dDenLdPn*NumL) / (DenL*DenL));
			A[0][2] = (lCoefficient.LINE_SCALE / lCoefficient.HEIGHT_SCALE)*((dNumLdHn*DenL - dDenLdHn*NumL) / (DenL*DenL));
			A[1][0] = (lCoefficient.SAMP_SCALE / lCoefficient.LONG_SCALE)*((dNumSdLn*DenS - dDenSdLn*NumS) / (DenS*DenS));
			A[1][1] = (lCoefficient.SAMP_SCALE / lCoefficient.LAT_SCALE)*((dNumSdPn*DenS - dDenSdPn*NumS) / (DenS*DenS));
			A[1][2] = (lCoefficient.SAMP_SCALE / lCoefficient.HEIGHT_SCALE)*((dNumSdHn*DenS - dDenSdHn*NumS) / (DenS*DenS));

			Ltmp = ((NumL / DenL)*lCoefficient.LINE_SCALE + lCoefficient.LINE_OFF);
			Stmp = ((NumS / DenS)*lCoefficient.SAMP_SCALE + lCoefficient.SAMP_OFF);

			double L = lAffine.lineb0 + lAffine.lineb1*Stmp + lAffine.lineb2*Ltmp;
			double S = lAffine.samplea0 + lAffine.samplea1*Stmp + lAffine.samplea2*Ltmp;

			l[0][0] = lpt[i].line - L;
			l[1][0] = lpt[i].sample - S;

			//……………………处理右片………………………………

			Pn = (ptObj[i].P - rCoefficient.LAT_OFF) / rCoefficient.LAT_SCALE;
			Ln = (ptObj[i].L - rCoefficient.LONG_OFF) / rCoefficient.LONG_SCALE;
			Hn = (ptObj[i].H - rCoefficient.HEIGHT_OFF) / rCoefficient.HEIGHT_SCALE;


			NumL = getaccumulation1(rCoefficient.LINE_NUM_COEFF, Ln, Pn, Hn);
			DenL = getaccumulation1(rCoefficient.LINE_DEN_COEFF, Ln, Pn, Hn);
			NumS = getaccumulation1(rCoefficient.SAMP_NUM_COEFF, Ln, Pn, Hn);
			DenS = getaccumulation1(rCoefficient.SAMP_DEN_COEFF, Ln, Pn, Hn);


			dNumLdLn = getpartialderivativeofL1(rCoefficient.LINE_NUM_COEFF, Ln, Pn, Hn);
			dNumLdPn = getpartialderivativeofP1(rCoefficient.LINE_NUM_COEFF, Ln, Pn, Hn);
			dNumLdHn = getpartialderivativeofH1(rCoefficient.LINE_NUM_COEFF, Ln, Pn, Hn);

			dDenLdLn = getpartialderivativeofL1(rCoefficient.LINE_DEN_COEFF, Ln, Pn, Hn);
			dDenLdPn = getpartialderivativeofP1(rCoefficient.LINE_DEN_COEFF, Ln, Pn, Hn);
			dDenLdHn = getpartialderivativeofH1(rCoefficient.LINE_DEN_COEFF, Ln, Pn, Hn);

			dNumSdLn = getpartialderivativeofL1(rCoefficient.SAMP_NUM_COEFF, Ln, Pn, Hn);
			dNumSdPn = getpartialderivativeofP1(rCoefficient.SAMP_NUM_COEFF, Ln, Pn, Hn);
			dNumSdHn = getpartialderivativeofH1(rCoefficient.SAMP_NUM_COEFF, Ln, Pn, Hn);

			dDenSdLn = getpartialderivativeofL1(rCoefficient.SAMP_DEN_COEFF, Ln, Pn, Hn);
			dDenSdPn = getpartialderivativeofP1(rCoefficient.SAMP_DEN_COEFF, Ln, Pn, Hn);
			dDenSdHn = getpartialderivativeofH1(rCoefficient.SAMP_DEN_COEFF, Ln, Pn, Hn);

			A[2][0] = (rCoefficient.LINE_SCALE / rCoefficient.LONG_SCALE)*((dNumLdLn*DenL - dDenLdLn*NumL) / (DenL*DenL));
			A[2][1] = (rCoefficient.LINE_SCALE / rCoefficient.LAT_SCALE)*((dNumLdPn*DenL - dDenLdPn*NumL) / (DenL*DenL));
			A[2][2] = (rCoefficient.LINE_SCALE / rCoefficient.HEIGHT_SCALE)*((dNumLdHn*DenL - dDenLdHn*NumL) / (DenL*DenL));
			A[3][0] = (rCoefficient.SAMP_SCALE / rCoefficient.LONG_SCALE)*((dNumSdLn*DenS - dDenSdLn*NumS) / (DenS*DenS));
			A[3][1] = (rCoefficient.SAMP_SCALE / rCoefficient.LAT_SCALE)*((dNumSdPn*DenS - dDenSdPn*NumS) / (DenS*DenS));
			A[3][2] = (rCoefficient.SAMP_SCALE / rCoefficient.HEIGHT_SCALE)*((dNumSdHn*DenS - dDenSdHn*NumS) / (DenS*DenS));


			Ltmp = ((NumL / DenL)*rCoefficient.LINE_SCALE + rCoefficient.LINE_OFF);
			Stmp = ((NumS / DenS)*rCoefficient.SAMP_SCALE + rCoefficient.SAMP_OFF);

			L = rAffine.lineb0 + rAffine.lineb1*Stmp + rAffine.lineb2*Ltmp;
			S = rAffine.samplea0 + rAffine.samplea1*Stmp + rAffine.samplea2*Ltmp;

			l[2][0] = rpt[i].line - L;
			l[3][0] = rpt[i].sample - S;


			for (int p = 0; p<4; p++)
			{
				for (int q = 0; q<3; q++)
				{
					AT[q][p] = A[p][q];
				}
			}

			matrixmulti(*AT, *A, *ATA, 3, 4, 3);
			matrixmulti(*AT, *l, *ATl, 3, 4, 1);

			if (inv(*ATA, 3))
			{
				matrixmulti(*ATA, *ATl, *d, 3, 3, 1);

			}

			else break;

			ptObj[i].L += d[0][0];

			ptObj[i].P += d[1][0];

			ptObj[i].H += d[2][0];

			numofiterative++;

		}//end of while

		 //lpoint<<i<<ptl.r[i]<<ptl.c[i]<<endl;
		numofiterative = 0;
		d[0][0] = 1;
		d[1][0] = 1;
		d[2][0] = 50;

	}//end of for
}


//求逆矩阵
int GetEpipolarImage::inv(double *m1, int n)
{
	int *is, *js;
	int i, j, k, l, u, v;
	double temp, max_v;
	is = (int *)malloc(n * sizeof(int));
	js = (int *)malloc(n * sizeof(int));
	if (is == NULL || js == NULL)
	{
		printf("out of memory!\n");
		return(0);
	}
	for (k = 0; k<n; k++)
	{
		max_v = 0.0;
		for (i = k; i<n; i++)
			for (j = k; j<n; j++)
			{
				temp = fabs(m1[i*n + j]);
				if (temp>max_v)
				{
					max_v = temp; is[k] = i; js[k] = j;
				}
			}
		if (max_v == 0.0)
		{
			free(is); free(js);
			printf("invers is not availble!\n");
			return(0);
		}
		if (is[k] != k)
			for (j = 0; j<n; j++)
			{
				u = k*n + j; v = is[k] * n + j;
				temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
			}
		if (js[k] != k)
			for (i = 0; i<n; i++)
			{
				u = i*n + k; v = i*n + js[k];
				temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
			}
		l = k*n + k;
		m1[l] = 1.0 / m1[l];
		for (j = 0; j<n; j++)
			if (j != k)
			{
				u = k*n + j;
				m1[u] *= m1[l];
			}
		for (i = 0; i<n; i++)
			if (i != k)
				for (j = 0; j<n; j++)
					if (j != k)
					{
						u = i*n + j;
						m1[u] -= m1[i*n + k] * m1[k*n + j];
					}
		for (i = 0; i<n; i++)
			if (i != k)
			{
				u = i*n + k;
				m1[u] *= -m1[l];
			}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j<n; j++)
			{
				u = k*n + j; v = js[k] * n + j;
				temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
			}
		if (is[k] != k)
			for (i = 0; i<n; i++)
			{
				u = i*n + k; v = i*n + is[k];
				temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
			}
	}
	free(is);
	free(js);
	return(1);
}

//乘
void GetEpipolarImage::matrixmulti(double *r1, double *r2, double *r, int m, int n, int p)
{
	for (int i = 0; i<m*p; i++)
		*(r + i) = 0;
	for (int i = 0; i<m; i++)
	{
		for (int j = 0; j<p; j++)
		{
			for (int k = 0; k<n; k++)
			{
				*(r + i*p + j) += *(r1 + i*n + k)*(*(r2 + k*p + j));
			}

		}
	}
}

double GetEpipolarImage::getaccumulation1(double *rpc, double L, double P, double H)
{
	double M[20];
	double S = 0;
	M[0] = 1, M[1] = L, M[2] = P, M[3] = H,
		M[4] = L*P, M[5] = L*H, M[6] = P*H, M[7] = L*L, M[8] = P*P, M[9] = H*H,
		M[10] = P*L*H, M[11] = L*L*L, M[12] = L*P*P, M[13] = L*H*H, M[14] = L*L*P,
		M[15] = P*P*P, M[16] = P*H*H, M[17] = L*L*H, M[18] = P*P*H, M[19] = H*H*H;

	for (int j = 0; j<20; j++)
	{
		S += rpc[j] * M[j];
	}
	return S;

}

double GetEpipolarImage::getaccumulation1(double rpc[20], SATPoint3D & ObjPt)
{
	double M[20];
	double S = 0;
	double L = ObjPt.L;
	double P = ObjPt.P;
	double H = ObjPt.H;
	M[0] = 1, M[1] = L, M[2] = P, M[3] = H,
		M[4] = L*P, M[5] = L*H, M[6] = P*H, M[7] = L*L, M[8] = P*P, M[9] = H*H,
		M[10] = P*L*H, M[11] = L*L*L, M[12] = L*P*P, M[13] = L*H*H, M[14] = L*L*P,
		M[15] = P*P*P, M[16] = P*H*H, M[17] = L*L*H, M[18] = P*P*H, M[19] = H*H*H;

	for (int j = 0; j<20; j++)
	{
		S += rpc[j] * M[j];
	}
	return S;

}
//获取仿射变换的逆矩阵
// 	L0=lineb0+lineb1*L+lineb2*S      //右边地面坐标的投影
//  S0=samplea0+samplea1*S+samplea2*L //原始影像上的坐标
//  L = dstlineb0+dstlineb1*L0+dstlineb2*S0
//  S = dstsamplea0+dstsamplea1*L0+dstsamplea2*s0
bool GetEpipolarImage::GetInverseAffPara(RPCImAffine &srcAffPara, RPCImAffine &dstAffPara)
{
	double b0 = srcAffPara.lineb0, b1 = srcAffPara.lineb1, b2 = srcAffPara.lineb2;
	double a0 = srcAffPara.samplea0, a1 = srcAffPara.samplea1, a2 = srcAffPara.samplea2;

	double fenmu = a1*b2 - a2*b1;

	if (fenmu == 0)
	{
		return false;
	}
	//
	//
	double inverA[6];
	//Line 相关
	inverA[0] = (a2*b0 - a0*b2) / fenmu;
	inverA[1] = b2 / fenmu;
	inverA[2] = a2 / fenmu;


	//SAMPLE 相关
	inverA[3] = (a0*b1 - a1*b0) / fenmu;
	inverA[4] = -b1 / fenmu;
	inverA[5] = a1 / fenmu;

	dstAffPara.samplea0 = inverA[0];
	dstAffPara.samplea1 = inverA[1];
	dstAffPara.samplea2 = inverA[2];

	dstAffPara.lineb0 = inverA[3];
	dstAffPara.lineb1 = inverA[4];
	dstAffPara.lineb2 = inverA[5];

	return true;
}

double GetEpipolarImage::getpartialderivativeofL1(double *rpc, double L, double P, double H)
{
	return(rpc[1] + rpc[4] * P + rpc[5] * H + 2 * rpc[7] * L + rpc[10] * P*H + 3 * rpc[11] * L*L
		+ rpc[12] * P*P + rpc[13] * H*H + 2 * rpc[14] * P*L + 2 * rpc[17] * L*H);
}
double GetEpipolarImage::getpartialderivativeofP1(double *rpc, double L, double P, double H)
{
	return(rpc[2] + rpc[4] * L + rpc[6] * H + 2 * rpc[8] * P + rpc[10] * L*H + 2 * rpc[12] * L*P
		+ rpc[14] * L*L + 3 * rpc[15] * P*P + rpc[16] * H*H + 2 * rpc[18] * P*H);
}
double GetEpipolarImage::getpartialderivativeofH1(double *rpc, double L, double P, double H)
{
	return(rpc[3] + rpc[5] * L + rpc[6] * P + 2 * rpc[9] * H + rpc[10] * P*L + 2 * rpc[13] * L*H
		+ 2 * rpc[16] * P*H + rpc[17] * L*L + rpc[18] * P*P + 3 * rpc[19] * H*H);
}

void GetEpipolarImage::calOriginalCoord1(Mat disp, vector<vector<float>>kb,Mat &x, Mat &y)
{
	//反算原图像坐标
	int sampleCount = disp.cols;
	int lineCount = disp.rows;
	x = Mat::zeros(sampleCount*lineCount, 1, CV_32F);
	y = Mat::zeros(sampleCount*lineCount, 1, CV_32F);
	for (int i = 0; i < sampleCount; i++)
	{
		for (int j = 0; j < lineCount; j++)
		{
			if (disp.at<float>(j, i) != 0)
			{
				int lx = lineCount - 1 - i, ly = j;
				int rx = lx - disp.at<float>(ly, lx), ry = ly;

				double centerX = sampleCount / 2.0;
				double centerY = lineCount / 2.0;

				int lx_90 = (ly - centerY) + centerX;
				int ly_90 = -(lx - centerX) + centerY;
				int rx_90 = (ry - centerY) + centerX;
				int ry_90 = -(rx - centerX) + centerY;

				CvPoint tmppoint;
				double lk = kb[lx_90][0];
				double lb = kb[lx_90][1];
				double rk = kb[lx_90][2];
				double rb = kb[lx_90][3];

				int originalLY = ly_90;
				int originalLX = ly_90 / lk - lb / lk;
				int originalRY = ry_90;
				int originalRX = ry_90 / rk - rb / rk;
				x.at<float>(i*lineCount + j) = originalRX - originalLX;
				y.at<float>(i*lineCount + j) = originalRY - originalLY;
			}
		}
	}
}
//创建核线影像
void GetEpipolarImage::createEpipolarImage1()
{
	IplImage *IplLeftEpipolarImage, *IplRightEpipolarImage;
	IplLeftEpipolarImage = cvCreateImage(cvSize(sampleCount, lineCount), IPL_DEPTH_8U, 1);
	IplRightEpipolarImage = cvCreateImage(cvSize(sampleCount, lineCount), IPL_DEPTH_8U, 1);

	uchar* leftData = (uchar*)IplLeftEpipolarImage->imageData;
	uchar* rightData = (uchar*)IplRightEpipolarImage->imageData;

	char fileName[] = "E:/linekb.txt";
	ofstream fout(fileName);

	for (int tmpSample = 1; tmpSample <= sampleCount; tmpSample++) {
		SATPoint2D tmpPoint;
		tmpPoint.sample = tmpSample;
		tmpPoint.line = lineCount / 2;
		calCounterpartLines(tmpPoint, 20);

		fout << tmpSample << " " << lk << " " << lb << " " << rk << " " << rb << endl;

		for (int tmpLine = 1; tmpLine <= lineCount; tmpLine++) {
			int lTmpX = tmpLine / lk - lb / lk;
			if (lTmpX <= sampleCount&&lTmpX > 0) {
				leftEpipolarImage[(tmpLine - 1)*sampleCount + (tmpSample - 1)] = leftImage[(tmpLine - 1)*sampleCount + (lTmpX - 1)];
				int test = leftImage[(tmpLine - 1)*sampleCount + (lTmpX - 1)];
				leftData[(tmpLine - 1)*sampleCount + (tmpSample - 1)] = leftImage[(tmpLine - 1)*sampleCount + (lTmpX - 1)];

			}

			int rTmpX = tmpLine / rk - rb / rk;
			if (rTmpX <= sampleCount&&rTmpX > 0) {
				rightEpipolarImage[(tmpLine - 1)*sampleCount + (tmpSample - 1)] = rightImage[(tmpLine - 1)*sampleCount + (rTmpX - 1)];
				int test = rightImage[(tmpLine - 1)*sampleCount + (rTmpX - 1)];
				rightData[(tmpLine - 1)*sampleCount + (tmpSample - 1)] = rightImage[(tmpLine - 1)*sampleCount + (rTmpX - 1)];
			}
		}
	}
	Mat le1 = cvarrToMat(IplLeftEpipolarImage);
	Mat re1 = cvarrToMat(IplRightEpipolarImage);
	cv::imwrite("leftimg2.tif", le1);
	cv::imwrite("rightimg2.tif", re1);
	cv::namedWindow("IplLeftEpipolarImage", 0);
	cv::namedWindow("IplRightEpipolarImage", 0);
	cv::imshow("IplLeftEpipolarImage", le1);
	cv::imshow("IplRightEpipolarImage", re1);

	//写文件
	cv::waitKey();
}
void GetEpipolarImage::setoffset(double x, double y)
{
	offx = x;
	offy = y;
}