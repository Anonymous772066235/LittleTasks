#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "BasicFunction.h"
#include "GetEpipolarImage.h"
#include <SemiGlobalMatching.h>

using namespace std;

class CSamePointTrackSGM
{
private:
	std::vector<std::string>  m_imglist;
	string m_outpath;
	float m_mindisp, m_maxdisp;
	string m_projpath;
	string m_logfilename;

	vector<vector<ImagePoint>> m_samepoints;

	int m_gridr, m_gridc;//用于追踪同名点时格网间隔

private:
	//读取RPC
	RPCcoeffcient readrpc(char* filename, cv::Mat rpcdata);
	bool isPointInRect(cv::Point p, cv::Rect rect);
	//生成格网
	void CMeshGrid(const cv::Range xr, const cv::Range yr, cv::Mat& outX, cv::Mat& outY);
	//分块的密集匹配
	cv::Mat GetDispSGM(cv::Mat le, cv::Mat re, float mindisp, float maxdisp);
	cv::Mat GetDispSGM(cv::Mat le, cv::Mat re, float mindisp, float maxdisp, int patchnum,int i_);
	//根据视差图追踪同名点,带有格网参数用于点抽稀
	vector<vector<ImagePoint>> SamePointTrack(vector<cv::Mat> dispxs, vector<cv::Mat> dispys, int gridr, int gridc);
	//视差图转位移
	void calOriginalCoord(cv::Mat disp, vector<vector<float>>kb, cv::Mat& x, cv::Mat& y, double offx, double offy);
	//同名点核线影像-->原始影像
	vector<vector<ImagePoint>> HXImage2Image(cv::Mat disp, vector<vector<float>>kb, double offx, double offy, int gridr, int gridc, int id);
	//根据范围筛选同名点
	vector<vector<ImagePoint>> SamePointChoose(vector<vector<ImagePoint>> ps, cv::Rect roi);

public:
	CSamePointTrackSGM();
	~CSamePointTrackSGM();
	CSamePointTrackSGM(std::string projpath,std::vector<std::string> imglistpath, string imgpath, string rpcpath,
		float mind, float maxd,
		 int gridr = 10, int gridc = 10)
	{
		create( projpath,imglistpath,  mind, maxd,  gridr, gridc);
	}
	//构造
	void create(std::string projpath,std::vector<std::string> imglistpath,
		float mind, float maxd,
		 int gridr = 10, int gridc = 10);
	//影像list路径
	void SetImageListPath(std::vector<std::string> imglistpath)
	{
		m_imglist= imglistpath;
	
	}


	//最小视差
	void SetMinDisp(float mind)
	{
		m_mindisp = mind;
	}
	float GetMinDisp()
	{
		return m_mindisp;
	}
	//最大视差
	void SetMaxDisp(float maxd)
	{
		m_maxdisp = maxd;
	}
	float GetMaxDisp()
	{
		return m_maxdisp;
	}
	//工程路径
	void SetProjectPath(string projpath)
	{
		m_projpath = projpath;
	}
	string GetProjectPath()
	{
		return m_projpath;
	}
	//采样间隔R
	void SetGridR(int gridr)
	{
		m_gridr = gridr;
	}
	int GetGridR()
	{
		return m_gridr;
	}
	//采样间隔R
	void SetGridC(int gridc)
	{
		m_gridc = gridc;
	}
	int GetGridC()
	{
		return m_gridc;
	}
	//计算
	void compute();
	//获取结果
	vector<vector<ImagePoint>> GetSamePoints()
	{
		return m_samepoints;
	}
};

