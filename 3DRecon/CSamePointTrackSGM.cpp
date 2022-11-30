#include "CSamePointTrackSGM.h"


CSamePointTrackSGM::CSamePointTrackSGM()
{
	m_imglist = vector<string> {""};
	m_mindisp = 0;
	m_maxdisp = 0;
	// m_projpath = " ";

	m_gridr = 10;
	m_gridc = 10;
}

CSamePointTrackSGM::~CSamePointTrackSGM()
{

}

void CSamePointTrackSGM::create(std::string projpath,std::vector<std::string> imglist, 
	float mind, float maxd,
	int gridr, int gridc)
{
	/// @brief 
	/// @param imglist 影像list路径
	/// @param mind 最小视差
	/// @param maxd 最大视差
	/// @param projpath 工程路径
	/// @param gridr 行格网大小
	/// @param gridc 列格网大小

	m_imglist = imglist;

	m_mindisp = mind;
	m_maxdisp = maxd;

	m_gridr = gridr;
	m_gridc = gridc;
	m_projpath=projpath;
}




//vector<cv::Point2f> ppppp1, ppppp2;



void CSamePointTrackSGM::compute()
{
	//开始流程化处理生成视差图
	//vector<cv::Mat> all_disps;
	vector<cv::Mat> dispxs(m_imglist.size() - 1), dispys(m_imglist.size() - 1);
	cv::Mat tmp;
	GetEpipolarImage gei;
	//vector<vector<vector<float>>> all_kbs;

	for (size_t i = 0; i < m_imglist.size() - 1; i++)
	{



		//读取影像-------------------------------------------------------------
		string imgpath1 =m_imglist[i] + ".tif";
		string imgpath2 =m_imglist[i + 1] + ".tif";
		cv::Mat img1 = BasicFunction::ReadTifoneBand(imgpath1);
		cv::Mat img2 = BasicFunction::ReadTifoneBand(imgpath2);
		/*cv::Mat img1 = cv::imread(imgpath1, 0);
		cv::Mat img2 = cv::imread(imgpath2, 0);
		img1 = BasicFunction::LinearStrech(img1, 0.01);
		img2 = BasicFunction::LinearStrech(img2, 0.01);*/
		//有些传感器有毛病，影像大小不一样,这里需要统一大小，把右下角裁掉
		int rows = img1.rows <= img2.rows ? img1.rows : img2.rows;
		int cols = img1.cols <= img2.cols ? img1.cols : img2.cols;
		img1 = img1(cv::Rect(0, 0, cols, rows));
		img2 = img2(cv::Rect(0, 0, cols, rows));


		/*cv::namedWindow("le", 0);
		cv::imshow("le", img1);
		cv::namedWindow("re", 0);
		cv::imshow("re", img2);
		cv::waitKey();*/

		//读取RPC----------------------------------------------------------------
		string rpcpath1, rpcpath2 ;
		
		if(access((m_imglist[i] + "_rpc.txt").c_str(),0)==0)
		{
			rpcpath1 =  m_imglist[i] + "_rpc.txt";
			rpcpath2 = m_imglist[i + 1] + "_rpc.txt";

		}else if(access((m_imglist[i] + "_RPC.txt").c_str(),0)==0){
			rpcpath1 = m_imglist[i] + "_RPC.txt";
			rpcpath2 = m_imglist[i + 1] + "_RPC.txt";
		}


		cout<<endl<<rpcpath1<<endl;
		cout<<endl<<rpcpath2<<endl;

		RPCcoeffcient rpc1 = readrpc((char*)rpcpath1.c_str(), tmp);
		RPCcoeffcient rpc2 = readrpc((char*)rpcpath2.c_str(), tmp);

		cout<<endl<<"Epipolar image start"<<endl;

		//生成核线影像-----------------------------------------------------------
		gei.setLRImage(img1, img2);
		gei.setRPCcoefAndRPCImAffine(rpc1, rpc2);
		gei.setoffset(0, 0);
		cv::Mat le, re;
		vector<vector<float>> kb;
		kb = gei.createEpipolarImage(le, re);
		//all_kbs.push_back(kb);
		cout<<endl<<"Epipolar image done"<<endl;

		//test
		/*cv::namedWindow("le", 0);
		cv::imshow("le", le);
		cv::namedWindow("re", 0);
		cv::imshow("re", re);
		cv::waitKey();*/
		//cv::imwrite("le_cg.jpg", le);
		//cv::imwrite("re_cg.jpg", re);

		//密集匹配----------------------------------------------------------------
		//这里转置是可选的.考虑到纠正出的核线是上下关系
		cv::transpose(le, le);
		cv::transpose(re, re);

		std::string RunningStateMsg=BasicFunction::Convert2Json_new(200,"3DReconstruction is running",0.8*double(i)/double(m_imglist.size()-1),"false");
		BasicFunction::Write2Text(m_projpath+"/RunningState_Log.txt",RunningStateMsg);

		cv::Mat disp_float = GetDispSGM(le, re, m_mindisp, m_maxdisp, 255,i);


		//玩砸了，翻车了，只能两两匹配了，不能追踪了，就这样吧---------------------------------------
		vector<vector<ImagePoint>> tps = HXImage2Image(-disp_float, kb, 0, 0, m_gridr, m_gridc, i);
		

		m_samepoints.insert(m_samepoints.end(), tps.begin(), tps.end());

		//test
		// int num(200);
		// cv::Rect rect(num, num, img1.cols - 2 * num, img1.rows - 2 * num);
		// for (int i(0); i < tps.size(); i += 100)
		// {
		// 	cv::Point2f pp1 = cv::Point2f(tps[i][0].X, tps[i][0].Y);
		// 	cv::Point2f pp2 = cv::Point2f(tps[i][1].X, tps[i][1].Y);
		// 	if (isPointInRect(pp1, rect) && isPointInRect(pp2, rect))
		// 	{
		// 		cv::Mat mimg1 = img1(cv::Rect(int(pp1.x - num / 2), int(pp1.y - num / 2), num, num)).clone();
		// 		cv::Mat mimg2 = img2(cv::Rect(int(pp2.x - num / 2), int(pp2.y - num / 2), num, num)).clone();
		// 		cv::Mat mimg;
		// 		cv::hconcat(mimg1, mimg2, mimg);
		// 		cv::imwrite("test_chang/"
		// 			+ to_string(i) + ".jpg", mimg);
		// 	}
		// }

		//test2
		//画
		// for (int i(0); i < m_samepoints.size(); ++i)
		// {
		// 	cv::circle(img1, cv::Point(m_samepoints[i][0].X, m_samepoints[i][0].Y), 2, 255);
		// 	cv::circle(img2, cv::Point(m_samepoints[i][1].X, m_samepoints[i][1].Y), 2, 255);
		// }
		// //存
		// cv::imwrite("imgl.png", img1);
		// cv::imwrite("imgr.png", img1);


		//test
		/*int num(200);
		cv::Rect rect(num, num, img1.cols - 2 * num, img1.rows - 2 * num);
		vector<cv::Mat> tdx, tdy;
		tdx.push_back(dispxs[i]);
		tdy.push_back(dispys[i]);*/
		/*for (int i(0); i < ppppp1.size(); i += 100)
		{
			cv::Point2f pp1 = ppppp1[i];
			cv::Point2f pp2 = ppppp2[i];
			if (isPointInRect(pp1, rect) && isPointInRect(pp2, rect))
			{
				cv::Mat mimg1 = img1(cv::Rect(int(pp1.x - num / 2), int(pp1.y - num / 2), num, num)).clone();
				cv::Mat mimg2 = img2(cv::Rect(int(pp2.x - num / 2), int(pp2.y - num / 2), num, num)).clone();
				cv::Mat mimg;
				cv::hconcat(mimg1, mimg2, mimg);
				cv::imwrite("E:/aidaza/yuanwang/SGMfengzahung/SamePointTrack/SamePointTrack/test_chang/"
					+ to_string(i) + ".jpg", mimg);
			}
		}*/
		/*vector<vector<ImagePoint>> tp = SamePointTrack(tdx, tdy, m_gridr, m_gridc);
		for (int i(0); i<tp.size(); ++i)
		{
			cv::Point2f pp1(tp[i][0].X, tp[i][0].Y);
			cv::Point2f pp2(tp[i][1].X, tp[i][1].Y);
			if (isPointInRect(pp1, rect) && isPointInRect(pp2, rect))
			{
				cv::Mat mimg1 = img1(cv::Rect(int(pp1.x - num / 2), int(pp1.y - num / 2), num, num)).clone();
				cv::Mat mimg2 = img2(cv::Rect(int(pp2.x - num / 2), int(pp2.y - num / 2), num, num)).clone();
				cv::Mat mimg;
				cv::hconcat(mimg1, mimg2, mimg);
				cv::imwrite("E:/aidaza/yuanwang/SGMfengzahung/SamePointTrack/SamePointTrack/test_chang/"
					+ to_string(i) + ".jpg", mimg);
			}
		}*/

		//入栈
		//all_disps.push_back(-disp_float);//位移是反过来的！！！


		//test-------
		/*int num(200);
		cv::Rect rect(num, num, le.cols - 2 * num, le.rows - 2 * num);
		vector<cv::Mat> td;
		td.push_back(-disp_float);
		vector<vector<ImagePoint>> tp = SamePointTrack(td, m_gridr, m_gridc);
		for (int i(0); i<tp.size(); ++i)
		{
			cv::Point2f pp1(tp[i][0].X, tp[i][0].Y);
			cv::Point2f pp2(tp[i][1].X, tp[i][1].Y);
			if (isPointInRect(pp1, rect) && isPointInRect(pp2, rect))
			{
				cv::Mat mimg1 = le(cv::Rect(int(pp1.x - num / 2), int(pp1.y - num / 2), num, num)).clone();
				cv::Mat mimg2 = re(cv::Rect(int(pp2.x - num / 2), int(pp2.y - num / 2), num, num)).clone();
				cv::Mat mimg;
				cv::hconcat(mimg1, mimg2, mimg);
				cv::imwrite("E:/aidaza/yuanwang/SGMfengzahung/SamePointTrack/SamePointTrack/test/"
					+ to_string(i) + ".jpg", mimg);
			}
		}*/

		//test2
		/*cv::transpose(le, le);
		cv::transpose(re, re);
		int num(200);
		cv::Rect rect(num, num, le.cols - 2 * num, le.rows - 2 * num);
		vector<cv::Mat> td;
		td.push_back(-disp_float);
		vector<vector<ImagePoint>> tp = SamePointTrack(td, m_gridr, m_gridc);
		for (int i(0); i < tp.size(); ++i)
		{
			cv::Point2f pp1(tp[i][0].Y, tp[i][0].X);
			cv::Point2f pp2(tp[i][1].Y, tp[i][1].X);
			if (isPointInRect(pp1, rect) && isPointInRect(pp2, rect))
			{
				cv::Mat mimg1 = le(cv::Rect(int(pp1.x - num / 2), int(pp1.y - num / 2), num, num)).clone();
				cv::Mat mimg2 = re(cv::Rect(int(pp2.x - num / 2), int(pp2.y - num / 2), num, num)).clone();
				cv::Mat mimg;
				cv::hconcat(mimg1, mimg2, mimg);
				cv::imwrite("E:/aidaza/yuanwang/SGMfengzahung/SamePointTrack/SamePointTrack/test_taishan/"
					+ to_string(i) + ".jpg", mimg);
			}
		}*/

	
	}

	/*cv::Rect rect = cv::Rect(0, 2200, 12000, 700);
	m_samepoints = SamePointChoose(m_samepoints, rect);*/

	//根据视差图进行同名点追踪------------------------------------------------------------------
	//m_samepoints = SamePointTrack(dispxs, dispys, m_gridr, m_gridc);
	//test
	//读取影像-------------------------------------------------------------
	/*vector<cv::Mat> timgs(m_imglist.size());
	for (size_t i = 0; i < m_imglist.size(); i++)
	{
		string imgpath = m_imgpath + "/" + m_imglist[i] + ".tif";
		timgs[i] = BasicFunction::ReadTifoneBand(imgpath);
	}*/
	//画
	/*for (int i(0); i < m_samepoints.size(); ++i)
	{
		for (int j(0); j < m_samepoints[i].size(); ++j)
		{
			cv::circle(timgs[m_samepoints[i][j].ID], cv::Point(m_samepoints[i][j].X, m_samepoints[i][j].Y), 2, 255);
		}
	}*/
	//存
	/*for (size_t i = 0; i < m_imglist.size(); i++)
	{
		cv::imwrite("test_plot/" + to_string(i) + ".png", timgs[i]);
	}*/




}


bool CSamePointTrackSGM::isPointInRect(cv::Point p, cv::Rect rect)
{
	int x0 = p.x;
	int y0 = p.y;
	int x1 = rect.x;
	int y1 = rect.y;
	int x2 = rect.x + rect.width;
	int y2 = rect.y + rect.height;
	if ((x1 < x0 && x0 < x2) &&
		(y1 < y0 && y0 < y2))
		return true;
	else
		return false;
}

cv::Mat CSamePointTrackSGM::GetDispSGM(cv::Mat le, cv::Mat re, float mindisp, float maxdisp)
{
	SemiGlobalMatching::SGMOption sgm_option;
	// 聚合路径数
	sgm_option.num_paths = 8;
	// 候选视差范围
	sgm_option.min_disparity = mindisp;
	sgm_option.max_disparity = maxdisp;
	// census窗口类型
	sgm_option.census_size = SemiGlobalMatching::Census5x5;
	// 一致性检查
	sgm_option.is_check_lr = true;
	sgm_option.lrcheck_thres = 1.0f;
	// 唯一性约束
	sgm_option.is_check_unique = true;
	sgm_option.uniqueness_ratio = 0.99;
	// 剔除小连通区
	sgm_option.is_remove_speckles = true;
	sgm_option.min_speckle_aera = 40;
	// 惩罚项P1、P2
	sgm_option.p1 = 10;
	sgm_option.p2_init = 150;
	// 视差图填充
	sgm_option.is_fill_holes = false;
	//匹配
	SemiGlobalMatching sgm;
	const sint32 width = static_cast<uint32>(le.cols);
	const sint32 height = static_cast<uint32>(le.rows);
	sgm.Initialize(width, height, sgm_option);
	float* disp = new float[uint32(width * height)];
	sgm.Match(le.data, re.data, disp);
	cv::Mat disp_mat(height, width, CV_32FC1, disp);
	//disp_mat.setTo(0, disp_mat == Invalid_Float);
	cv::Mat disp_return = disp_mat.clone();
	delete[]disp;
	return disp_return;
}

cv::Mat CSamePointTrackSGM::GetDispSGM(cv::Mat le, cv::Mat re, float mindisp, float maxdisp, int patchnum,int i_)
{
	//核线是这样的
	/*--------------------        --------------------
	   --------------------        --------------------
		--------------------        --------------------
		 --------------------        --------------------*/
		 //首先按行分块-------------------------------------
	int rows = le.rows;
	int cols = le.cols;
	int patchlen = rows / patchnum;
	cout << "\nrow:\t"<<rows <<endl;
	cout << "col:\t"<<cols<<endl<<endl;;

	std::vector<cv::Mat> limgs, rimgs;
	for (size_t i = 0; i < rows; i += patchlen)
	{
		if(i+patchlen<rows)
		{
			limgs.push_back(le(cv::Rect(0, i, cols, patchlen)).clone());
			rimgs.push_back(re(cv::Rect(0, i, cols, patchlen)).clone());
		}
		else
		{
			limgs.push_back(le(cv::Rect(0, i, cols, rows-i)).clone());
			rimgs.push_back(re(cv::Rect(0, i, cols, rows-i)).clone());
		}
		
	}
	//然后进行密集匹配---------------------------------
	//我想多线程，但是狗玩意SGM太耗内存了，所以慢慢来吧
	cv::Mat all_disp = cv::Mat::zeros(rows, cols, CV_32FC1);
	for (size_t i = 0; i < limgs.size(); i++)
	// for (size_t i = 0; i < 1; i++)
	{
		cout << "\nThe " << to_string(i + 1) << " of " << to_string(limgs.size()) << " blocks" << endl;

		std::string RunningStateMsg=BasicFunction::Convert2Json_new(200,"3DReconstruction is running",0.8*double(i_)/double(m_imglist.size()-1)+(double(i+1)/double(limgs.size()))*(0.8/double(m_imglist.size()-1)),"false");
		BasicFunction::Write2Text(m_projpath+"/RunningState_Log.txt",RunningStateMsg);

		cv::Mat disp_mat = GetDispSGM(limgs[i], rimgs[i], mindisp, maxdisp);

		if((i+1)*patchlen<rows)
		{
			disp_mat.copyTo(all_disp(cv::Rect(0, i * patchlen, cols, patchlen)));
		}
		else
		{
			disp_mat.copyTo(all_disp(cv::Rect(0, i * patchlen, cols, rows-i * patchlen)));
		}
	}

	//test
	// cv::Mat tdisp = all_disp.clone();
	// tdisp.setTo(0, all_disp == Invalid_Float);
	// cv::Mat disp_show = cv::Mat(tdisp.rows, tdisp.cols, CV_8UC1);
	// cv::normalize(tdisp, disp_show, 0, 255, cv::NORM_MINMAX, CV_8UC1);
	// cv::namedWindow("视差图", 0);
	// cv::imshow("视差图", disp_show);
	
	// cv::imwrite("视差图.jpg", disp_show);
	// cv::waitKey();

	return all_disp;
}

RPCcoeffcient CSamePointTrackSGM::readrpc(char* filename, cv::Mat rpcdata)
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

	cv::Mat temp(90, 1, CV_32F);

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

void CSamePointTrackSGM::CMeshGrid(const cv::Range xr, const cv::Range yr, cv::Mat& outX, cv::Mat& outY)
{
	vector<float> x, y;
	for (int i = xr.start; i <= xr.end; i++)
		x.push_back(i);
	for (int i = yr.start; i <= yr.end; i++)
		y.push_back(i);
	//需要对x进行一次转置
	repeat(cv::Mat(x).t(), y.size(), 1, outX);
	repeat(cv::Mat(y), 1, x.size(), outY);
}

vector<vector<ImagePoint>>  CSamePointTrackSGM::SamePointTrack(vector<cv::Mat> dispxs, vector<cv::Mat> dispys, int gridr, int gridc)
{
	vector<vector<ImagePoint>> samepoints;
	int rows = dispxs[0].rows;
	int cols = dispxs[0].cols;
	int imgnum = dispxs.size() + 1;
	//根据像素位移计算每个点在每一帧的位置-----------------
	//生成格网
	cv::Mat fx0, fy0;//首帧点
	CMeshGrid(cv::Range(0, cols - 1), cv::Range(0, rows - 1), fx0, fy0);
	//位移累加获取同名点位置
	vector<cv::Mat> fx, fy;
	fx.push_back(fx0);
	fy.push_back(fy0);
	cv::Mat tfx0 = fx0.clone();
	cv::Mat tfy0 = fy0.clone();
	for (int i(0); i < imgnum - 1; ++i)
	{
		//位移累加
		cv::Mat tfx, tfy;
		tfx.create(rows, cols, CV_32FC1);
		tfy.create(rows, cols, CV_32FC1);
		add(tfx0, dispxs[i], tfx);
		add(tfy0, dispys[i], tfy);
		tfy = tfy0.clone();
		//累加后位移入栈
		fx.push_back(tfx);
		fy.push_back(tfy);
		//更新首帧
		tfx0 = tfx.clone();
		tfy0 = tfy.clone();
	}
	//同名点追踪---------------------------------------------------
	vector<ImagePoint> pbunchs;//同名点串
	//第一帧开始每个点追踪
	float x(0), y(0);
	//int lajiit(1);
	for (int i(0); i < rows; i += gridr)
	{
		for (int j(0); j < cols; j += gridc)
		{
			pbunchs.clear();
			//每帧循环进行点追踪------------------
			for (int k(0); k < imgnum; ++k)
			{
				x = fx[k].at<float>(i, j);
				y = fy[k].at<float>(i, j);
				if (i == 2499 && j == 2499)
				{
					cout << x << "  " << y << endl;
				}
				//仅保存影像范围内、且牛皮的同名点
				if ((0 <= x && x < cols) &&
					(0 <= y && y < rows) &&
					x != Invalid_Float && y != Invalid_Float)
				{
					pbunchs.push_back(ImagePoint(x, y, k));
				}
			}
			//对于追踪到的同名点进行保存--------
			if (pbunchs.size() > 1)
			{
				samepoints.push_back(pbunchs);
			}
		}
	}
	return samepoints;
}

void CSamePointTrackSGM::calOriginalCoord(cv::Mat disp, vector<vector<float>>kb, cv::Mat& x, cv::Mat& y, double offx, double offy)
{
	//反算原图像坐标
	int sampleCount = disp.cols;
	int lineCount = disp.rows;
	x = Invalid_Float * cv::Mat::ones(sampleCount, lineCount, CV_32FC1);
	y = Invalid_Float * cv::Mat::ones(sampleCount, lineCount, CV_32FC1);
	for (int i = 0; i < lineCount - 1; i++)
	{
		for (int j = 0; j < sampleCount; j++)
		{
			if (disp.at<float>(i, j) != Invalid_Float)
			{
				float lx = j;
				float ly = i;
				float rx = j + disp.at<float>(i, j);
				float ry = i;

				//反算
				double lk = kb[i][0];
				double lb = kb[i][1];
				double rk = kb[i][2];
				double rb = kb[i][3];

				double originalLY = lx;
				double originalLX = (lx + offy) / lk - lb / lk - offx;
				double originalRY = rx;
				double originalRX = (rx + offy) / rk - rb / rk - offx;

				/*if (i == 2499 && j == 2499)
				{
					originalRY = 2499;
					originalRX = (2499 + offy) / rk - rb / rk - offx;
					cout << originalLY << "  " << originalLX << endl;
					cout << originalRY << "  " << originalRX << endl;
				}*/

				if (originalLX > 0 && originalLX < lineCount && originalLY>0 && originalLY < sampleCount)
				{
					x.at<float>(int(originalLY + .5), int(originalLX + .5)) = originalRX - originalLX;//按列进行存储，可根据需要进行调整
					y.at<float>(int(originalLY + .5), int(originalLX + .5)) = originalRY - originalLY;//按列进行存储，可根据需要进行调整

					//ppppp1.push_back(cv::Point2f(originalLX, originalLY));
					//ppppp2.push_back(cv::Point2f(originalRX, originalRY));
				}
			}
		}
	}
}



vector<vector<ImagePoint>> CSamePointTrackSGM::HXImage2Image(cv::Mat disp, vector<vector<float>>kb, double offx, double offy, int gridr, int gridc, int id)
{
	//反算原图像坐标
	int sampleCount = disp.cols;
	int lineCount = disp.rows;
	vector<vector<ImagePoint>> all_ps;
	vector<ImagePoint> tps;
	for (int i = 0; i < lineCount - 1; i += gridr)
	{
		for (int j = 0; j < sampleCount; j += gridc)
		{
			tps.clear();
			if (!isinf(disp.at<float>(i, j)))
			{
				float lx = j;
				float ly = i;
				float rx = j + disp.at<float>(i, j);
				float ry = i;

				//反算
				double lk = kb[i][0];
				double lb = kb[i][1];
				double rk = kb[i][2];
				double rb = kb[i][3];

				double originalLY = lx;
				double originalLX = (lx + offy) / lk - lb / lk - offx;
				double originalRY = rx;
				double originalRX = (rx + offy) / rk - rb / rk - offx;

				if (originalLX > 0 && originalLX < lineCount && originalLY>0 && originalLY < sampleCount)
				{
					tps.push_back(ImagePoint(originalLX, originalLY, id));
					tps.push_back(ImagePoint(originalRX, originalRY, id + 1));
					all_ps.push_back(tps);
				}
			}
		}
	}
	return all_ps;
}



vector<vector<ImagePoint>> CSamePointTrackSGM::SamePointChoose(vector<vector<ImagePoint>> ps, cv::Rect roi)
{
	vector<vector<ImagePoint>> nps;
	for (int i(0); i < ps.size(); ++i)
	{
		cv::Point p = cv::Point(ps[i][0].X, ps[i][0].Y);
		if (isPointInRect(p, roi))
		{
			nps.push_back(ps[i]);
		}
	}
	return nps;
}