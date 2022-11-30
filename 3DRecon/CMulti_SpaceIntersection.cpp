#include "CMulti_SpaceIntersection.h"
#include "COORCONV.h"
std::mutex mtx;

//构造、析构函数
CMulti_SpaceIntersection::CMulti_SpaceIntersection()
{

}
CMulti_SpaceIntersection::~CMulti_SpaceIntersection()
{

}
//设置其他参数
void CMulti_SpaceIntersection::SetOtherParams(int T_ProjectErr, int k, double m,std::string projpath)
{
	p_T_ProjectErr = T_ProjectErr;
	p_k = p_k;
	p_m = p_m;
	p_projpath=projpath;
}

//设置RPC参数
void CMulti_SpaceIntersection::SetInputPathParams( std::vector<std::string> ImageListPath)
{
	p_ImageList = ImageListPath;
}

//输入同名点
void CMulti_SpaceIntersection::SetInputCorreImagePoints(std::vector<std::vector<ImagePoint>> CorrespondingImagePoints)
{
	p_CorrespondingImagePoints = CorrespondingImagePoints;
}

//设置左右影像的RPC参数
void CMulti_SpaceIntersection::SetRPCParams(RPCcoeffcient& lRPCcoef, RPCcoeffcient& rRPCcoef)
{
	p_lRPC = lRPCcoef;
	p_rRPC = rRPCcoef;
}
void CMulti_SpaceIntersection::SetImgAffine(RPCImAffine lAffine, RPCImAffine rAffine)
{
	p_lAffine = lAffine;
	p_rAffine = rAffine;
}
//分割函数
void CMulti_SpaceIntersection::Stringsplit(const std::string& str, const std::string& splits, std::vector<std::string>& res)
{
	if (str == "")		return;
	//在字符串末尾也加入分隔符，方便截取最后一段
	std::string strs = str + splits;
	size_t pos = strs.find(splits);
	int step = splits.size();

	// 若找不到内容则字符串搜索函数返回 npos
	while (pos != strs.npos)
	{
		std::string temp = strs.substr(0, pos);
		res.push_back(temp);
		//去掉已分割的字符串,在剩下的字符串中进行分割
		strs = strs.substr(pos + step, strs.size());
		pos = strs.find(splits);
	}
}
std::vector<std::vector<double>> CMulti_SpaceIntersection::ReadRPCImAffineParamsFromTxt(std::string FilePath)
{
	std::ifstream ifile(FilePath);//读取文件
	std::streampos len = ifile.tellg();//获取文件长度

	std::vector<std::vector<double>> data;
	if (ifile) //文件打开是否成功
	{
		std::string str; //temp
		std::cout << "file open scessful" << std::endl;

		while (std::getline(ifile, str))
		{
			std::vector<std::string> strtmp;
			std::vector<double> datatmp;
			Stringsplit(str, "	    ", strtmp);
			for (int i = 0; i < strtmp.size(); i++)
			{
				datatmp.push_back(std::stof(strtmp[i]));
			}
			data.push_back(datatmp);
		}
		ifile.close();
		return data;
	}
	else
	{
		std::vector<std::vector<double>> RPCAdj(2,std::vector<double>(3,0));
		RPCAdj[0][1]=1;
		RPCAdj[1][2]=1;
		return RPCAdj;

	}
}
RPCcoeffcient CMulti_SpaceIntersection::ReadRPC(char* filename)
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

	for (int i = 0; i < 20; i++)
	{
		rpcdata1.LINE_DEN_COEFF[i] = D_L[i];
		rpcdata1.LINE_NUM_COEFF[i] = N_L[i];
		rpcdata1.SAMP_DEN_COEFF[i] = D_S[i];
		rpcdata1.SAMP_NUM_COEFF[i] = N_S[i];
	}
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

//RPC参数及其平差后的参数
bool CMulti_SpaceIntersection::LoadRPCParams()
{
	if (p_ImageList.size() != 0)
	{
		for (int i = 0; i < p_ImageList.size(); i++)
		{
			//std::vector<std::vector<double>> RPCAdj = ReadRPCImAffineParamsFromTxt(p_RPCImAffineParamsPaths[i]);
			std::vector<std::vector<double>> RPCAdj = ReadRPCImAffineParamsFromTxt((char*)( p_ImageList[i] + "_rpc-r.adj").c_str());

			std::cout<<"RPCAdj"<<std::endl;

			std::cout<<RPCAdj[0][0]<<" "<<RPCAdj[0][1]<<" "<<RPCAdj[0][2]<<std::endl;
			std::cout<<RPCAdj[1][0]<<" "<<RPCAdj[1][1]<<" "<<RPCAdj[1][2]<<std::endl;

		

			if (RPCAdj.size() != 0)
			{
				RPCImAffine iRPCImAffineParams;
				iRPCImAffineParams.lineb0 = RPCAdj[1][0];
				iRPCImAffineParams.lineb1 = RPCAdj[1][1];
				iRPCImAffineParams.lineb2 = RPCAdj[1][2];
				iRPCImAffineParams.samplea0 = RPCAdj[0][0];
				iRPCImAffineParams.samplea1 = RPCAdj[0][1];
				iRPCImAffineParams.samplea2 = RPCAdj[0][2];
				p_ImAffineParams.push_back(iRPCImAffineParams);
			}
			else
			{
				RPCImAffine iRPCImAffineParams;
				iRPCImAffineParams.lineb0 = RPCAdj[1][0];
				iRPCImAffineParams.lineb1 = RPCAdj[1][1];
				iRPCImAffineParams.lineb2 = RPCAdj[1][2];
				iRPCImAffineParams.samplea0 = RPCAdj[0][0];
				iRPCImAffineParams.samplea1 = RPCAdj[0][1];
				iRPCImAffineParams.samplea2 = RPCAdj[0][2];
				p_ImAffineParams.push_back(iRPCImAffineParams);
			}
			//RPCcoeffcient iRPCParams = ReadRPC((char *)p_RPCParamsPaths[i].data());

			string rpcpath;
			if(access((p_ImageList[i]+ "_rpc.txt").c_str(),0)==0)
			{
				rpcpath =  p_ImageList[i] + "_rpc.txt";

			}else if(access((p_ImageList[i] + "_RPC.txt").c_str(),0)==0){
				rpcpath =  p_ImageList[i] + "_RPC.txt";
			}

			RPCcoeffcient iRPCParams = ReadRPC((char*)(rpcpath.data()));
			p_RPCParams.push_back(iRPCParams);
		}
		return true;
	}
	else
	{
		return false;
	}
}
//保存文件
namespace
{
	template<typename _T>
	static void writeStream(std::ofstream& file, const _T* data, int rows, int cols, int matStep, int cn)
	{
		int x, y;
		for (y = 0; y < rows; y++)
		{
			const _T* pData = data + y * matStep;
			for (x = 0; x < cols * cn; x += cn)
			{
				if (cn == 1)
				{
					file << pData[x] << " ";//单通道数据，每列数据用 隔开
					continue;
				}
				for (int channel = 0; channel < cn; channel++)
					file << pData[x + channel] << (channel + 1 < cn ? "," : "\t");//多通道数据，同一列不同通道用','隔开，									//每列数据用Tab隔开
			}
			file << std::endl;
		}
	}
}


void saveAsText(std::string filename, cv::InputArray _src)
{

	std::ofstream outFile(filename.c_str(), std::ios_base::out);
	if (!outFile.is_open())
	{
		std::cout << "Fail to open file (" << filename << ")" << std::endl;
		return;
	}
	cv::Mat src = _src.getMat();//此方法不支持vector<Mat>或vector<vector<>>，具体可参考OpenCV源码(core/src/matrix.cpp)
	
	CV_Assert(src.empty() == false);

	int depth = src.depth();
	int matStep = (int)(src.step / src.elemSize1());
	const uchar* data = src.data;
	int cn = src.channels();
	
	int rows = src.rows;
	int cols = src.cols;

	if (depth == CV_8U)
		writeStream(outFile, (const uchar*)data, rows, cols, matStep, cn);
	else if (depth == CV_8S)
		writeStream(outFile, (const schar*)data, rows, cols, matStep, cn);
	else if (depth == CV_16U)
		writeStream(outFile, (const ushort*)data, rows, cols, matStep, cn);
	else if (depth == CV_16S)
		writeStream(outFile, (const short*)data, rows, cols, matStep, cn);
	else if (depth == CV_32S)
		writeStream(outFile, (const int*)data, rows, cols, matStep, cn);
	else if (depth == CV_32F)
	{
		std::streamsize pp = outFile.precision();
		outFile.precision(16);//控制输出精度，可以自行修改
		writeStream(outFile, (const float*)data, rows, cols, matStep, cn);
		outFile.precision(pp);//恢复默认输出精度
	}
	else if (depth == CV_64F)
	{
		std::streamsize pp = outFile.precision();
		outFile.precision(16);
		writeStream(outFile, (const double*)data, rows, cols, matStep, cn);
		outFile.precision(pp);
	}
	else	
	{

	}
}

//坐标转换
cv::Point3f CMulti_SpaceIntersection::BLH2XYZ(double b, double l, double h)
{
	//const double alpha = 1 / 298.257223563;
	//const double a = 6378137;
	//double b1, l1;//转化度分秒为度和弧度

	//b1 = b / 180 * acos(-1);//得到弧度
	//l1 = l / 180 * acos(-1);//得到弧度
	//double e = sqrt(2 * alpha - alpha * alpha);
	//double N = a / sqrt(1 - e * e * sin(b1) * sin(b1));
	//cv::Point3f xyz;
	//xyz.x = (N + h) * cos(b1) * cos(l1);
	//xyz.y = (N + h) * cos(b1) * sin(l1);
	//xyz.z = (N * (1 - e * e) + h) * sin(b1);
	cv::Point3f xyz;
	double x, y;
	LonLat2UTM(l, b, x, y);
	xyz.x = x;
	xyz.y = y;
	xyz.z = h;
	return xyz;
}
SATPoint3D CMulti_SpaceIntersection::XYZ2BLH(double x, double y, double z)
{
	const double alpha = 1 / 298.257223563;
	const double a = 6378137;
	double b, h, N, H = 0, pcs, b1, l1, l2;//设置b、l、h代替BLH，求取精度
	double eps = 0.00001;//精度为0.00001
	b = 0; h = 0;
	double e = sqrt(2 * alpha - alpha * alpha);
	b = atan(z / (sqrt(x * x + y * y) * (1 - e * e)));
	l1 = atan2(y, x);//利用atan2函数，得到弧度
					 //根据精度需要求取b和h
	do
	{
		N = a / sqrt(1 - e * e * sin(b) * sin(b));
		b = atan(z / (sqrt(x * x + y * y) * (1 - e * e * N / (N + h))));
		h = sqrt(x * x + y * y) / cos(b) - N;
		pcs = h - H;
		H = h;
	} while (fabs(pcs) >= eps);
	//将弧度转化为度分秒形式
	b1 = b / acos(-1) * 180;
	l2 = l1 / acos(-1) * 180;
	SATPoint3D blh;
	blh.P = b1;
	blh.L = l2;
	blh.H = H;
	return blh;
}

vector<cv::Point3f> CMulti_SpaceIntersection::SOR(vector<cv::Point3f> InputPointCloud, int kdnum, double m)
{

	//8线程
	int T_PtNum = InputPointCloud.size() / 8;
	// int T_PtNum = InputPointCloud.size() / 8;

	// //第1个线程
	vector<cv::Point3f> v_pt3d1(InputPointCloud.begin(), InputPointCloud.begin() + T_PtNum);
	vector<double> dis1;
	future<void> ft1 = async(std::launch::async, [&]
		{
			cv::Mat m_pt3d = cv::Mat(v_pt3d1).reshape(1);
			//trees参数设置为2
			cv::flann::KDTreeIndexParams indexParams(2);
			 //kd树索引建
			cv::flann::Index kdtree;
			kdtree.build(m_pt3d, cv::flann::KDTreeIndexParams(1), cvflann::FLANN_DIST_EUCLIDEAN);

			vector<float> vecQuery(3);
			vector<int> vecIndex(kdnum + 1);
			vector<float> vecDist(kdnum + 1);
			
			//搜索参数
			cv::flann::SearchParams params(8);

			for (int i = 0; i < v_pt3d1.size(); i++)
			{

				vecQuery[0] = v_pt3d1[i].x;
				vecQuery[1] = v_pt3d1[i].y;
				vecQuery[2] = v_pt3d1[i].z;
				kdtree.knnSearch(vecQuery, vecIndex, vecDist, kdnum + 1, params);
				double distance_sum = 0;
				for (int j = 1; j < kdnum + 1; j++)
				{
					distance_sum = distance_sum + sqrt(vecDist[j]);
				}
				double average = distance_sum / kdnum;
				dis1.push_back(average);
			}
		});

	//第2个线程
	vector<cv::Point3f> v_pt3d2(InputPointCloud.begin() + T_PtNum, InputPointCloud.begin() + 2 * T_PtNum);
	vector<double> dis2;
	future<void> ft2 = async(std::launch::async, [&]
		{
			cv::Mat m_pt3d = cv::Mat(v_pt3d2).reshape(1);
			//trees参数设置为2
			cv::flann::KDTreeIndexParams indexParams(2);
			//kd树索引建
			cv::flann::Index kdtree;
			kdtree.build(m_pt3d, cv::flann::KDTreeIndexParams(1), cvflann::FLANN_DIST_EUCLIDEAN);

			vector<float> vecQuery(3);
			vector<int> vecIndex(kdnum + 1);
			vector<float> vecDist(kdnum + 1);

			//搜索参数
			cv::flann::SearchParams params(8);

			for (int i = 0; i < v_pt3d2.size(); i++)
			{

				vecQuery[0] = v_pt3d2[i].x;
				vecQuery[1] = v_pt3d2[i].y;
				vecQuery[2] = v_pt3d2[i].z;
				kdtree.knnSearch(vecQuery, vecIndex, vecDist, kdnum + 1, params);
				double distance_sum = 0;
				for (int j = 1; j < kdnum + 1; j++)
				{
					distance_sum = distance_sum + sqrt(vecDist[j]);
				}
				double average = distance_sum / kdnum;
				dis2.push_back(average);
			}
		});

	//第3个线程
	vector<cv::Point3f> v_pt3d3(InputPointCloud.begin() + 2 * T_PtNum, InputPointCloud.begin() + 3 * T_PtNum);
	vector<double> dis3;
	future<void> ft3 = async(std::launch::async, [&]
		{
			cv::Mat m_pt3d = cv::Mat(v_pt3d3).reshape(1);
			//trees参数设置为2
			cv::flann::KDTreeIndexParams indexParams(2);
			//kd树索引建
			cv::flann::Index kdtree;
			kdtree.build(m_pt3d, cv::flann::KDTreeIndexParams(1), cvflann::FLANN_DIST_EUCLIDEAN);

			vector<float> vecQuery(3);
			vector<int> vecIndex(kdnum + 1);
			vector<float> vecDist(kdnum + 1);

			//搜索参数
			cv::flann::SearchParams params(8);

			for (int i = 0; i < v_pt3d3.size(); i++)
			{

				vecQuery[0] = v_pt3d3[i].x;
				vecQuery[1] = v_pt3d3[i].y;
				vecQuery[2] = v_pt3d3[i].z;
				kdtree.knnSearch(vecQuery, vecIndex, vecDist, kdnum + 1, params);
				double distance_sum = 0;
				for (int j = 1; j < kdnum + 1; j++)
				{
					distance_sum = distance_sum + sqrt(vecDist[j]);
				}
				double average = distance_sum / kdnum;
				dis3.push_back(average);
			}
		});
	//第4个线程
	vector<cv::Point3f> v_pt3d4(InputPointCloud.begin()+ 3* T_PtNum, InputPointCloud.begin() + 4*T_PtNum);
	vector<double> dis4;
	future<void> ft4 = async(std::launch::async, [&]
		{
			cv::Mat m_pt3d = cv::Mat(v_pt3d4).reshape(1);
			//trees参数设置为2
			cv::flann::KDTreeIndexParams indexParams(2);
			//kd树索引建
			cv::flann::Index kdtree;
			kdtree.build(m_pt3d, cv::flann::KDTreeIndexParams(1), cvflann::FLANN_DIST_EUCLIDEAN);

			vector<float> vecQuery(3);
			vector<int> vecIndex(kdnum + 1);
			vector<float> vecDist(kdnum + 1);

			//搜索参数
			cv::flann::SearchParams params(8);

			for (int i = 0; i < v_pt3d4.size(); i++)
			{

				vecQuery[0] = v_pt3d4[i].x;
				vecQuery[1] = v_pt3d4[i].y;
				vecQuery[2] = v_pt3d4[i].z;
				kdtree.knnSearch(vecQuery, vecIndex, vecDist, kdnum + 1, params);
				double distance_sum = 0;
				for (int j = 1; j < kdnum + 1; j++)
				{
					distance_sum = distance_sum + sqrt(vecDist[j]);
				}
				double average = distance_sum / kdnum;
				dis4.push_back(average);
			}
		});
	//第5个线程
	vector<cv::Point3f> v_pt3d5(InputPointCloud.begin() + 4 * T_PtNum, InputPointCloud.begin() + 5*T_PtNum);
	vector<double> dis5;
	future<void> ft5 = async(std::launch::async, [&]
		{
			cv::Mat m_pt3d = cv::Mat(v_pt3d5).reshape(1);
			//trees参数设置为2
			cv::flann::KDTreeIndexParams indexParams(2);
			//kd树索引建
			cv::flann::Index kdtree;
			kdtree.build(m_pt3d, cv::flann::KDTreeIndexParams(1), cvflann::FLANN_DIST_EUCLIDEAN);

			vector<float> vecQuery(3);
			vector<int> vecIndex(kdnum + 1);
			vector<float> vecDist(kdnum + 1);

			//搜索参数
			cv::flann::SearchParams params(8);

			for (int i = 0; i < v_pt3d5.size(); i++)
			{

				vecQuery[0] = v_pt3d5[i].x;
				vecQuery[1] = v_pt3d5[i].y;
				vecQuery[2] = v_pt3d5[i].z;
				kdtree.knnSearch(vecQuery, vecIndex, vecDist, kdnum + 1, params);
				double distance_sum = 0;
				for (int j = 1; j < kdnum + 1; j++)
				{
					distance_sum = distance_sum + sqrt(vecDist[j]);
				}
				double average = distance_sum / kdnum;
				dis5.push_back(average);
			}
		});
	//第6个线程
	vector<cv::Point3f> v_pt3d6(InputPointCloud.begin() + 5 * T_PtNum, InputPointCloud.begin() +6* T_PtNum);
	vector<double> dis6;
	future<void> ft6 = async(std::launch::async, [&]
		{
			cv::Mat m_pt3d = cv::Mat(v_pt3d6).reshape(1);
			//trees参数设置为2
			cv::flann::KDTreeIndexParams indexParams(2);
			//kd树索引建
			cv::flann::Index kdtree;
			kdtree.build(m_pt3d, cv::flann::KDTreeIndexParams(1), cvflann::FLANN_DIST_EUCLIDEAN);

			vector<float> vecQuery(3);
			vector<int> vecIndex(kdnum + 1);
			vector<float> vecDist(kdnum + 1);

			//搜索参数
			cv::flann::SearchParams params(8);

			for (int i = 0; i < v_pt3d6.size(); i++)
			{

				vecQuery[0] = v_pt3d6[i].x;
				vecQuery[1] = v_pt3d6[i].y;
				vecQuery[2] = v_pt3d6[i].z;
				kdtree.knnSearch(vecQuery, vecIndex, vecDist, kdnum + 1, params);
				double distance_sum = 0;
				for (int j = 1; j < kdnum + 1; j++)
				{
					distance_sum = distance_sum + sqrt(vecDist[j]);
				}
				double average = distance_sum / kdnum;
				dis6.push_back(average);
			}
		});
	//第7个线程
	vector<cv::Point3f> v_pt3d7(InputPointCloud.begin() + 7 * T_PtNum, InputPointCloud.begin() + 8*T_PtNum);
	vector<double> dis7;
	future<void> ft7 = async(std::launch::async, [&]
		{
			cv::Mat m_pt3d = cv::Mat(v_pt3d7).reshape(1);
			//trees参数设置为2
			cv::flann::KDTreeIndexParams indexParams(2);
			//kd树索引建
			cv::flann::Index kdtree;
			kdtree.build(m_pt3d, cv::flann::KDTreeIndexParams(1), cvflann::FLANN_DIST_EUCLIDEAN);

			vector<float> vecQuery(3);
			vector<int> vecIndex(kdnum + 1);
			vector<float> vecDist(kdnum + 1);

			//搜索参数
			cv::flann::SearchParams params(8);

			for (int i = 0; i < v_pt3d7.size(); i++)
			{

				vecQuery[0] = v_pt3d7[i].x;
				vecQuery[1] = v_pt3d7[i].y;
				vecQuery[2] = v_pt3d7[i].z;
				kdtree.knnSearch(vecQuery, vecIndex, vecDist, kdnum + 1, params);
				double distance_sum = 0;
				for (int j = 1; j < kdnum + 1; j++)
				{
					distance_sum = distance_sum + sqrt(vecDist[j]);
				}
				double average = distance_sum / kdnum;
				dis7.push_back(average);
			}
		});
	//第8个线程
	vector<cv::Point3f> v_pt3d8(InputPointCloud.begin() + 8 * T_PtNum, InputPointCloud.end());
	vector<double> dis8;
	future<void> ft8 = async(std::launch::async, [&]
		{
			cv::Mat m_pt3d = cv::Mat(v_pt3d8).reshape(1);
			//trees参数设置为2
			cv::flann::KDTreeIndexParams indexParams(2);
			//kd树索引建
			cv::flann::Index kdtree;
			kdtree.build(m_pt3d, cv::flann::KDTreeIndexParams(1), cvflann::FLANN_DIST_EUCLIDEAN);

			vector<float> vecQuery(3);
			vector<int> vecIndex(kdnum + 1);
			vector<float> vecDist(kdnum + 1);

			//搜索参数
			cv::flann::SearchParams params(8);

			for (int i = 0; i < v_pt3d8.size(); i++)
			{

				vecQuery[0] = v_pt3d8[i].x;
				vecQuery[1] = v_pt3d8[i].y;
				vecQuery[2] = v_pt3d8[i].z;
				kdtree.knnSearch(vecQuery, vecIndex, vecDist, kdnum + 1, params);
				double distance_sum = 0;
				for (int j = 1; j < kdnum + 1; j++)
				{
					distance_sum = distance_sum + sqrt(vecDist[j]);
				}
				double average = distance_sum / kdnum;
				dis8.push_back(average);
			}
		});
	//线程等待
	ft1.wait();
	ft2.wait();
	ft3.wait();
	ft4.wait();
	ft5.wait();
	ft6.wait();
	ft7.wait();
	ft8.wait();
	// //pinjie
	vector<vector<double>>dis;
	dis.push_back(dis1);
	dis.push_back(dis2);
	dis.push_back(dis3);
	dis.push_back(dis4);
	dis.push_back(dis5);
	dis.push_back(dis6);
	dis.push_back(dis7);
	dis.push_back(dis8);


	


	//计算平均距离
	cv::Mat distance(InputPointCloud.size(), 1, CV_32F);
	int inum = 0;
	for (int i = 0; i < dis.size(); i++)
	{
		for (int j = 0; j < dis[i].size(); j++)
		{
			distance.at<float>(inum, 0) = dis[i][j];
			inum++;
		}
	}
	cv::Mat mean_distance;
	cv::Mat _std;


	cv::meanStdDev(distance, mean_distance, _std);

	//输出点云
	vector<cv::Point3f> OutputPointCloud;
	
	double t1 = mean_distance.at<double>(0, 0);
	double t2 = _std.at<double>(0, 0);
	double T = mean_distance.at<double>(0, 0) + m * _std.at<double>(0, 0);

	for (int i = 0; i < InputPointCloud.size(); i++)
	{
		if (distance.at<float>(i, 0) < T)
		{
			OutputPointCloud.push_back(InputPointCloud[i]);

		}
	}

	return OutputPointCloud;
}




//多片前方交会
bool  CMulti_SpaceIntersection::Multi_SpaceIntersection(std::vector<SATPoint3D>& PointCloud,std::string PointCloudSavePath)
{

	clock_t start, finish;
	double duration;
	using namespace std::chrono;

	start = clock();
	system_clock::time_point now = system_clock::now();
	std::time_t last = system_clock::to_time_t(now - std::chrono::hours(0));
	cout << endl<<std::put_time(std::localtime(&last), "%F %T") << " " << ":" << " " << "Start point cloud generation" << endl;
	//加载一下参数
	if (!LoadRPCParams())
	{
		return false;
	}



	// 在这里确定一下第一张影像的角点
	string xmlpath=p_ImageList[0]+".xml";
	vector<vector<double>> CornerPoints;
	BasicFunction::readXML4Boundary(xmlpath,CornerPoints);

	//输出为txt,另外算出平均经度以确定投影带号（默认为49）
	ofstream CorPfile;
	CorPfile.precision(16);
	CorPfile.setf(ios_base::showpoint);
	CorPfile.open(p_projpath+"/CornerPoints.txt");
	
	double average_Lon=0;
	int zone=49;
	if(CornerPoints.size()==4)
	{
		CorPfile  << "[";
		for (int i = 0; i <CornerPoints.size(); i++)
		{
			average_Lon+=CornerPoints[i][0]/4.0;
			CorPfile  << "[\""<<CornerPoints[i][0] << "\",\""<< CornerPoints[i][1]<<"\"]";
			if(i<CornerPoints.size()-1)
			{
				CorPfile  << ",";
			}
		}
		CorPfile  << "]";
		CorPfile.close();
		int zone=floor(average_Lon/6)+31;
		
	}
	else
	{
		cout << '\n' << "Something wrong in reading xml.\n";
	}


    std::cout<<"\naverage_Lon\t"<<average_Lon<<"\n";
	std::cout<<"\nzone\t"<<zone<<"\n";
	std::cout<<"\np_CorrespondingImagePoints.size()\t"<<p_CorrespondingImagePoints.size()<<"\n";

	
	int T_PtNum = p_CorrespondingImagePoints.size() / 8;

	//第1个线程
	//三维点纵坐标
	
	// vector<SATPoint3D>Points_BLH1;
	vector<cv::Point3f> Points_XYZ1;
	//Flag
	cv::Mat Flag1(T_PtNum, 1, CV_32S);
	//同名点
	vector<vector<ImagePoint>> CorImagePoints1(p_CorrespondingImagePoints.begin(), p_CorrespondingImagePoints.begin() + T_PtNum);
	//线程内
	future<void> ft1 = async(std::launch::async, [&]
		{
			UTMCoor xy;
			for (int i = 0; i < CorImagePoints1.size(); i++)
			{
				int iflag;
				SATPoint3D itmp_p3d;
				itmp_p3d = Multi_RPCIntersection(CorImagePoints1[i],iflag);
				// Points_BLH1.push_back(itmp_p3d);
				Flag1.at<int>(i, 0) = iflag;
				LatLonToUTMXY(DegToRad(itmp_p3d.P), DegToRad( itmp_p3d.L),zone,xy);
				cv::Point3f  ipt3d(xy.x,xy.y,itmp_p3d.H);
				Points_XYZ1.push_back(ipt3d);
			}
		});
	//第2个线程
	//三维点纵坐标
	// vector<SATPoint3D>Points_BLH2;
	vector<cv::Point3f> Points_XYZ2;
	//Flag
	cv::Mat Flag2(T_PtNum, 1, CV_32S);
	//同名点
	vector<vector<ImagePoint>> CorImagePoints2(p_CorrespondingImagePoints.begin()+ T_PtNum * 1, p_CorrespondingImagePoints.begin() + T_PtNum * 2);
	//线程内
	future<void> ft2 = async(std::launch::async, [&]
		{
			UTMCoor xy;
			for (int i = 0; i < CorImagePoints2.size(); i++)
			{
				int iflag;
				SATPoint3D itmp_p3d;
				itmp_p3d = Multi_RPCIntersection(CorImagePoints2[i], iflag);

				// Points_BLH2.push_back(itmp_p3d);
				Flag2.at<int>(i, 0) = iflag;
				LatLonToUTMXY(DegToRad(itmp_p3d.P), DegToRad( itmp_p3d.L),zone,xy);
				cv::Point3f  ipt3d(xy.x,xy.y,itmp_p3d.H);
				Points_XYZ2.push_back(ipt3d);
			}
		});
	//第3个线程
	//三维点纵坐标
	// vector<SATPoint3D>Points_BLH3;
	vector<cv::Point3f> Points_XYZ3;

	//Flag
	cv::Mat Flag3(T_PtNum, 1, CV_32S);
	//同名点
	vector<vector<ImagePoint>> CorImagePoints3(p_CorrespondingImagePoints.begin() + T_PtNum * 2, p_CorrespondingImagePoints.begin() + T_PtNum * 3);
	//线程内

	future<void> ft3 = async(std::launch::async, [&]
		{
			UTMCoor xy;
			for (int i = 0; i < CorImagePoints3.size(); i++)
			{
				int iflag;
				SATPoint3D itmp_p3d;
				itmp_p3d = Multi_RPCIntersection(CorImagePoints3[i], iflag);
				// Points_BLH3.push_back(itmp_p3d);
				Flag3.at<int>(i, 0) = iflag;
				LatLonToUTMXY(DegToRad(itmp_p3d.P), DegToRad( itmp_p3d.L),zone,xy);
				cv::Point3f  ipt3d(xy.x,xy.y,itmp_p3d.H);
				Points_XYZ3.push_back(ipt3d);
			}
		});
	//第4个线程
	//三维点纵坐标
	// vector<SATPoint3D>Points_BLH4;
	vector<cv::Point3f> Points_XYZ4;

	//Flag
	cv::Mat Flag4(T_PtNum, 1, CV_32S);
	//同名点
	vector<vector<ImagePoint>> CorImagePoints4(p_CorrespondingImagePoints.begin() + T_PtNum * 3, p_CorrespondingImagePoints.begin() + T_PtNum * 4);
	//线程内

	future<void> ft4 = async(std::launch::async, [&]
		{
			UTMCoor xy;
			for (int i = 0; i < CorImagePoints4.size(); i++)
			{
				int iflag;
				SATPoint3D itmp_p3d;
				itmp_p3d = Multi_RPCIntersection(CorImagePoints4[i], iflag);
				// Points_BLH4.push_back(itmp_p3d);
				Flag4.at<int>(i, 0) = iflag;
				LatLonToUTMXY(DegToRad(itmp_p3d.P), DegToRad( itmp_p3d.L),zone,xy);
				cv::Point3f  ipt3d(xy.x,xy.y,itmp_p3d.H);
				Points_XYZ4.push_back(ipt3d);
			}
		});
	//第5个线程
	//三维点纵坐标
	// vector<SATPoint3D>Points_BLH5;
	vector<cv::Point3f> Points_XYZ5;

	//Flag
	cv::Mat Flag5(T_PtNum, 1, CV_32S);
	//同名点
	vector<vector<ImagePoint>> CorImagePoints5(p_CorrespondingImagePoints.begin() + T_PtNum * 4, p_CorrespondingImagePoints.begin() + T_PtNum * 5);
	//线程内

	future<void> ft5 = async(std::launch::async, [&]
		{
			UTMCoor xy;
			for (int i = 0; i < CorImagePoints5.size(); i++)
			{
				int iflag;
				SATPoint3D itmp_p3d;
				itmp_p3d = Multi_RPCIntersection(CorImagePoints5[i], iflag);
				// Points_BLH5.push_back(itmp_p3d);
				Flag5.at<int>(i, 0) = iflag;
				LatLonToUTMXY(DegToRad(itmp_p3d.P), DegToRad( itmp_p3d.L),zone,xy);
				cv::Point3f  ipt3d(xy.x,xy.y,itmp_p3d.H);	
				Points_XYZ5.push_back(ipt3d);
			}
		});
	//第6个线程
	//三维点纵坐标
	// vector<SATPoint3D>Points_BLH6;
	vector<cv::Point3f> Points_XYZ6;

	//Flag
	cv::Mat Flag6(T_PtNum, 1, CV_32S);
	//同名点
	vector<vector<ImagePoint>> CorImagePoints6(p_CorrespondingImagePoints.begin() + T_PtNum * 5, p_CorrespondingImagePoints.begin() + T_PtNum * 6);
	//线程内

	future<void> ft6 = async(std::launch::async, [&]
		{
			UTMCoor xy;
			for (int i = 0; i < CorImagePoints6.size(); i++)
			{

				int iflag;
				SATPoint3D itmp_p3d;
				itmp_p3d = Multi_RPCIntersection(CorImagePoints6[i], iflag);
				// Points_BLH6.push_back(itmp_p3d);
				Flag6.at<int>(i, 0) = iflag;
				LatLonToUTMXY(DegToRad(itmp_p3d.P), DegToRad( itmp_p3d.L),zone,xy);
				cv::Point3f  ipt3d(xy.x,xy.y,itmp_p3d.H);	
				Points_XYZ6.push_back(ipt3d);
			}
		});
	//第7个线程
	//三维点纵坐标
	// vector<SATPoint3D>Points_BLH7;
	vector<cv::Point3f> Points_XYZ7;

	//Flag
	cv::Mat Flag7(T_PtNum, 1, CV_32S);
	//同名点
	vector<vector<ImagePoint>> CorImagePoints7(p_CorrespondingImagePoints.begin() + T_PtNum * 6, p_CorrespondingImagePoints.begin() + T_PtNum * 7);
	//线程内

	future<void> ft7 = async(std::launch::async, [&]
		{
			UTMCoor xy;
			for (int i = 0; i < CorImagePoints7.size(); i++)
			{
				int iflag;
				SATPoint3D itmp_p3d;
				itmp_p3d = Multi_RPCIntersection(CorImagePoints7[i], iflag);
				// Points_BLH7.push_back(itmp_p3d);
				Flag7.at<int>(i, 0) = iflag;
				LatLonToUTMXY(DegToRad(itmp_p3d.P), DegToRad( itmp_p3d.L),zone,xy);
				cv::Point3f  ipt3d(xy.x,xy.y,itmp_p3d.H);	
				Points_XYZ7.push_back(ipt3d);
			}
		});
	//第8个线程
	//三维点纵坐标
	// vector<SATPoint3D>Points_BLH8;
	vector<cv::Point3f> Points_XYZ8;

	//Flag
	cv::Mat Flag8(p_CorrespondingImagePoints.size() - T_PtNum*7, 1, CV_32S);
	//同名点
	vector<vector<ImagePoint>> CorImagePoints8(p_CorrespondingImagePoints.begin() + T_PtNum * 7, p_CorrespondingImagePoints.end());
	//线程内

	future<void> ft8 = async(std::launch::async, [&]
		{
			UTMCoor xy;
			for (int i = 0; i < CorImagePoints8.size(); i++)
			{

				int iflag;
				SATPoint3D itmp_p3d;
				itmp_p3d = Multi_RPCIntersection(CorImagePoints8[i], iflag);
				// Points_BLH8.push_back(itmp_p3d);
				Flag8.at<int>(i, 0) = iflag;
				LatLonToUTMXY(DegToRad(itmp_p3d.P), DegToRad( itmp_p3d.L),zone,xy);
				cv::Point3f  ipt3d(xy.x,xy.y,itmp_p3d.H);	
				Points_XYZ8.push_back(ipt3d);
			}
		});



	//线程等待
	ft1.wait();
	ft2.wait();
	ft3.wait();
	ft4.wait();
	ft5.wait();
	ft6.wait();
	ft7.wait();
	ft8.wait();


	cv::Mat Flag = cv::Mat::zeros(p_CorrespondingImagePoints.size(), 1, CV_32S);
	
	//Flag 拼接
	Flag1.copyTo(Flag.rowRange(0 * T_PtNum, 1 * T_PtNum));
	Flag2.copyTo(Flag.rowRange(1 * T_PtNum, 2 * T_PtNum));
	Flag3.copyTo(Flag.rowRange(2 * T_PtNum, 3 * T_PtNum));
	Flag4.copyTo(Flag.rowRange(3 * T_PtNum, 4 * T_PtNum));
	Flag5.copyTo(Flag.rowRange(4 * T_PtNum, 5 * T_PtNum));
	Flag6.copyTo(Flag.rowRange(5 * T_PtNum, 6 * T_PtNum));
	Flag7.copyTo(Flag.rowRange(6 * T_PtNum, 7 * T_PtNum));
	Flag8.copyTo(Flag.rowRange(7 * T_PtNum, Flag.rows));


	//三维点拼接
	std::vector<cv::Point3f> Points_XYZ;
	// std::vector<SATPoint3D> Points_BLH;


	for (int i = 0; i < T_PtNum; i++)
	{
		
		std::string RunningStateMsg=BasicFunction::Convert2Json_new(200,"3DReconstruction is running",0.8+0.1*(double(i+1)/double(T_PtNum)),"false");
		BasicFunction::Write2Text(p_projpath+"/RunningState_Log.txt",RunningStateMsg);
		//XYZ
		//前面七个线程的内容
		if (Flag1.at<int>(i, 0) == 1)
		{
			Points_XYZ.push_back(Points_XYZ1[i]);
			// Points_BLH.push_back(Points_BLH1[i]);
		}
		if (Flag2.at<int>(i, 0) == 1)
		{
			Points_XYZ.push_back(Points_XYZ2[i]);
			// Points_BLH.push_back(Points_BLH2[i]);
		}
		if (Flag3.at<int>(i, 0) == 1)
		{
			Points_XYZ.push_back(Points_XYZ3[i]);
			// Points_BLH.push_back(Points_BLH3[i]);
		}
		if (Flag4.at<int>(i, 0) == 1)
		{
			Points_XYZ.push_back(Points_XYZ4[i]);
			// Points_BLH.push_back(Points_BLH4[i]);
		}
		if (Flag5.at<int>(i, 0) == 1)
		{
			Points_XYZ.push_back(Points_XYZ5[i]);
			// Points_BLH.push_back(Points_BLH5[i]);
		}
		if (Flag6.at<int>(i, 0) == 1)
		{
			// Points_BLH.push_back(Points_BLH6[i]);
			Points_XYZ.push_back(Points_XYZ6[i]);
		}
		if (Flag7.at<int>(i, 0) == 1)
		{
			// Points_BLH.push_back(Points_BLH7[i]);
			Points_XYZ.push_back(Points_XYZ7[i]);
		}
		//最后一块
		if (i < Points_XYZ8.size())
		{

			if (Flag8.at<int>(i, 0) == 1)
			{
				Points_XYZ.push_back(Points_XYZ8[i]);
				// Points_BLH.push_back(Points_BLH8[i]);
			}
		}
	}


	// 极端情况下，出现SOR后点云为0的情况
	std::cout<<"Before SOR, point cloud size:"<<std::endl;
	std::cout<<Points_XYZ.size()<<std::endl;

	vector<cv::Point3f> Points_XYZ_filter = SOR(Points_XYZ, 50, 0.3);

	std::cout<<"After SOR, point cloud size:"<<std::endl;
	std::cout<<Points_XYZ_filter.size()<<std::endl;



	ofstream outputfile;
	outputfile.precision(16);
	outputfile.setf(ios_base::showpoint);
	outputfile.open(PointCloudSavePath);
	// 若滤波后的点的数量为0.则输出滤波前的点云
	if(Points_XYZ_filter.size()>0)
	{
		for (int i = 0; i <Points_XYZ_filter.size(); i++)
		{
			outputfile  << Points_XYZ_filter[i].x << " "<< Points_XYZ_filter[i].y<<" "<<Points_XYZ_filter[i].z<<endl;
		}
		outputfile.close();
	}
	else if(Points_XYZ_filter.size()==0)
	{
		cout << '\n' << "The filtered point could is null, output the orginal point cloud\n";
		for (int i = 0; i <Points_XYZ.size(); i++)
		{
			outputfile  << Points_XYZ[i].x << " "<< Points_XYZ[i].y<<" "<<Points_XYZ[i].z<<endl;
		}
		outputfile.close();

	}
	

	cout << '\n' << std::put_time(std::localtime(&last), "%F %T") << " " << ":" << " " << "Static target reconstruction completed";

}



cv::Mat CMulti_SpaceIntersection::GetAandRCandL(SATPoint2D* lpt, SATPoint3D* ptObj, RPCImAffine imaffine)
{
	//同名点数量
	int npt = 1;
	//flag是上面的循环次，将第一次和后面的次分开

	double A[2][3], l[2][1];
	/*d[0][0] = 1;
	d[1][0] = 1;
	d[2][0] = 50;*/
	double Pn = 0, Ln = 0, Hn = 0;
	double NumL = 0, DenL = 0, NumS = 0, DenS = 0;
	double dNumLdLn = 0.0, dNumLdPn = 0.0, dNumLdHn = 0.0;
	double dDenLdLn = 0.0, dDenLdPn = 0.0, dDenLdHn = 0.0;
	double dNumSdLn = 0.0, dNumSdPn = 0.0, dNumSdHn = 0.0;
	double dDenSdLn = 0.0, dDenSdPn = 0.0, dDenSdHn = 0.0;  //上面四行用来记录正解法偏导数

	cv::Mat ARCL(2 * npt, 5, CV_32F);
	double Ltmp, Stmp;

	for (int i = 0; i < npt; i++)
	{
		// while(numofiterative<30&&(fabs(d[0][0])>0.00001||fabs(d[1][0])>0.00001||fabs(d[2][0])>0.001))
		//……………………处理左片………………	
		Pn = (ptObj[i].P - p_lRPC.LAT_OFF) / p_lRPC.LAT_SCALE;
		Ln = (ptObj[i].L - p_lRPC.LONG_OFF) / p_lRPC.LONG_SCALE;
		Hn = (ptObj[i].H - p_lRPC.HEIGHT_OFF) / p_lRPC.HEIGHT_SCALE;

		NumL = getaccumulation1(p_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
		DenL = getaccumulation1(p_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
		NumS = getaccumulation1(p_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
		DenS = getaccumulation1(p_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);

		dNumLdLn = getpartialderivativeofL1(p_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
		dNumLdPn = getpartialderivativeofP1(p_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
		dNumLdHn = getpartialderivativeofH1(p_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);

		dDenLdLn = getpartialderivativeofL1(p_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
		dDenLdPn = getpartialderivativeofP1(p_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
		dDenLdHn = getpartialderivativeofH1(p_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);

		dNumSdLn = getpartialderivativeofL1(p_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
		dNumSdPn = getpartialderivativeofP1(p_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
		dNumSdHn = getpartialderivativeofH1(p_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);

		dDenSdLn = getpartialderivativeofL1(p_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);
		dDenSdPn = getpartialderivativeofP1(p_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);
		dDenSdHn = getpartialderivativeofH1(p_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);

		A[0][0] = (p_lRPC.LINE_SCALE / p_lRPC.LONG_SCALE) * ((dNumLdLn * DenL - dDenLdLn * NumL) / (DenL * DenL));
		A[0][1] = (p_lRPC.LINE_SCALE / p_lRPC.LAT_SCALE) * ((dNumLdPn * DenL - dDenLdPn * NumL) / (DenL * DenL));
		A[0][2] = (p_lRPC.LINE_SCALE / p_lRPC.HEIGHT_SCALE) * ((dNumLdHn * DenL - dDenLdHn * NumL) / (DenL * DenL));
		A[1][0] = (p_lRPC.SAMP_SCALE / p_lRPC.LONG_SCALE) * ((dNumSdLn * DenS - dDenSdLn * NumS) / (DenS * DenS));
		A[1][1] = (p_lRPC.SAMP_SCALE / p_lRPC.LAT_SCALE) * ((dNumSdPn * DenS - dDenSdPn * NumS) / (DenS * DenS));
		A[1][2] = (p_lRPC.SAMP_SCALE / p_lRPC.HEIGHT_SCALE) * ((dNumSdHn * DenS - dDenSdHn * NumS) / (DenS * DenS));

		Ltmp = ((NumL / DenL) * p_lRPC.LINE_SCALE + p_lRPC.LINE_OFF);
		Stmp = ((NumS / DenS) * p_lRPC.SAMP_SCALE + p_lRPC.SAMP_OFF);

		double L = imaffine.lineb0 + imaffine.lineb1 * Stmp + imaffine.lineb2 * Ltmp;
		double S = imaffine.samplea0 + imaffine.samplea1 * Stmp + imaffine.samplea2 * Ltmp;

		l[0][0] = lpt[i].line - L;
		l[1][0] = lpt[i].sample - S;

		//……………………处理右片………………………………

		ARCL.at<float>(0, 0) = A[0][0];
		ARCL.at<float>(0, 1) = A[0][1];
		ARCL.at<float>(0, 2) = A[0][2];
		ARCL.at<float>(1, 0) = A[1][0];
		ARCL.at<float>(1, 1) = A[1][1];
		ARCL.at<float>(1, 2) = A[1][2];

		ARCL.at<float>(0, 3) = Ltmp;
		ARCL.at<float>(1, 3) = Stmp;

		ARCL.at<float>(0, 4) = l[0][0];
		ARCL.at<float>(1, 4) = l[1][0];
		return ARCL;
	}
}

void CMulti_SpaceIntersection::LonLat2UTM(double longitude, double latitude, double& UTME, double& UTMN)
{
	double lat = latitude;
	double lon = longitude;

	double kD2R = pi / 180.0;
	double ZoneNumber = floor((lon - 1.5) / 3.0) + 1;
	double L0 = ZoneNumber * 3.0;

	double a = 6378137.0;
	double F = 298.257223563;
	double f = 1 / F;
	double b = a * (1 - f);
	double ee = (a * a - b * b) / (a * a);
	double e2 = (a * a - b * b) / (b * b);
	double n = (a - b) / (a + b);
	double n2 = (n * n);
	double n3 = (n2 * n);
	double n4 = (n2 * n2);
	double n5 = (n4 * n);
	double al = (a + b) * (1 + n2 / 4 + n4 / 64) / 2.0;
	double bt = -3 * n / 2 + 9 * n3 / 16 - 3 * n5 / 32.0;
	double gm = 15 * n2 / 16 - 15 * n4 / 32;
	double dt = -35 * n3 / 48 + 105 * n5 / 256;
	double ep = 315 * n4 / 512;
	double B = lat * kD2R;
	double L = lon * kD2R;
	L0 = L0 * kD2R;
	double l = L - L0;
	double cl = (cos(B) * l);
	double cl2 = (cl * cl);
	double cl3 = (cl2 * cl);
	double cl4 = (cl2 * cl2);
	double cl5 = (cl4 * cl);
	double cl6 = (cl5 * cl);
	double cl7 = (cl6 * cl);
	double cl8 = (cl4 * cl4);
	double lB = al * (B + bt * sin(2 * B) + gm * sin(4 * B) + dt * sin(6 * B) + ep * sin(8 * B));
	double t = tan(B);
	double t2 = (t * t);
	double t4 = (t2 * t2);
	double t6 = (t4 * t2);
	double Nn = a / sqrt(1 - ee * sin(B) * sin(B));
	double yt = e2 * cos(B) * cos(B);
	double N = lB;
	N = N + t * Nn * cl2 / 2;
	N = N + t * Nn * cl4 * (5 - t2 + 9 * yt + 4 * yt * yt) / 24;
	N = N + t * Nn * cl6 * (61 - 58 * t2 + t4 + 270 * yt - 330 * t2 * yt) / 720;
	N = N + t * Nn * cl8 * (1385 - 3111 * t2 + 543 * t4 - t6) / 40320;
	double E = Nn * cl;
	E = E + Nn * cl3 * (1 - t2 + yt) / 6;
	E = E + Nn * cl5 * (5 - 18 * t2 + t4 + 14 * yt - 58 * t2 * yt) / 120;
	E = E + Nn * cl7 * (61 - 479 * t2 + 179 * t4 - t6) / 5040;
	E = E + 500000;
	N = 0.9996 * N;
	E = 0.9996 * (E - 500000.0) + 500000.0;

	UTME = E;
	UTMN = N;
}

SATPoint3D CMulti_SpaceIntersection::Multi_RPCIntersection(vector<ImagePoint> sp, int& flag)
{
	mtx.lock();
	//先进行双像前方交会*********************************************************************
	int pnum = 1;
	int col = sp.size();
	SATPoint2D* lpt = new SATPoint2D[pnum];
	SATPoint2D* rpt = new SATPoint2D[pnum];
	SATPoint3D* ptObj = new SATPoint3D[pnum];

	//左片右片像点
	lpt[0].line = sp[0].Y;
	lpt[0].sample = sp[0].X;


	rpt[0].line = sp[col - 1].Y;
	rpt[0].sample = sp[col - 1].X;
	//使用第一张和最后一张	
	SetRPCParams(p_RPCParams[sp[0].ID], p_RPCParams[sp[col - 1].ID]);
	SetImgAffine(p_ImAffineParams[sp[0].ID], p_ImAffineParams[sp[col - 1].ID]);

	double lflag;

	if (RPCIntersection(lpt, rpt, pnum, ptObj))
	{

		lflag = 1;
	}
	else
	{
		lflag = 0;
	}
	//**************************************************************************************
	int n = 0;//迭代次数计数器
	while (n < 10)
	{
		cv::Mat AL = Multi_GetAandL(sp, ptObj);
		cv::Mat A(AL.rows, AL.cols - 1, CV_32F);
		cv::Mat L(AL.rows, 1, CV_32F);

		for (int i = 0; i < A.rows; i++)
		{
			for (int j = 0; j < A.cols; j++)
			{
				A.at<float>(i, j) = AL.at<float>(i, j);
			}
			L.at<float>(i, 0) = AL.at<float>(i, 3);
		}


		cv::Mat x;
		x = (A.t() * A).inv() * (A.t() * L);

		if (n == 9)
		{
			vector<float> vv;
			for (int j = 0; j < L.rows; j++)
			{
				vv.push_back(L.at<float>(j, 0));
			}
			for (int j = 0; j < L.rows; j++)
			{
				if (L.at<float>(j, 0) > 30 || ptObj[0].H < 0)
				{
					lflag = 0;
					break;
				}
			}

		}
		ptObj[0].L = ptObj[0].L + x.at<float>(0, 0);
		ptObj[0].P = ptObj[0].P + x.at<float>(1, 0);
		ptObj[0].H = ptObj[0].H + x.at<float>(2, 0);
		n = n + 1;

	}
	SATPoint3D pt3d;
	pt3d.P = ptObj[0].P;
	pt3d.L = ptObj[0].L;
	pt3d.H = ptObj[0].H;
	flag = lflag;
	mtx.unlock();
	return pt3d;

}
//返回多张影像的同名点的AL矩阵 flag1表达是否第一次迭代
cv::Mat CMulti_SpaceIntersection::Multi_GetAandL(vector<ImagePoint> pt2d, SATPoint3D* pt3d)
{

	//帧数
	int fnum = pt2d.size();

	//同名点数，基本为1

	int pnum = 1;

	//三维坐标初值已知

	SATPoint2D* lpt = new SATPoint2D[pnum];
	SATPoint2D* rpt = new SATPoint2D[pnum];

	cv::Mat AL(2 * fnum, 4, CV_32F);
	int jsq = 0;
	for (int i = 0; i < fnum; i++)
	{
		lpt[0].line = pt2d[i].Y;
		lpt[0].sample = pt2d[i].X;

		//设置RPC参数
		if (i <= fnum - 2)
		{
			SetRPCParams(p_RPCParams[pt2d[i].ID], p_RPCParams[pt2d[i].ID + 1]);
			SetImgAffine(p_ImAffineParams[pt2d[i].ID], p_ImAffineParams[pt2d[i].ID + 1]);
		}
		if (i == fnum - 1)
		{
			SetRPCParams(p_RPCParams[pt2d[i].ID], p_RPCParams[0]);
			SetImgAffine(p_ImAffineParams[pt2d[i].ID], p_ImAffineParams[pt2d[i].ID + 1]);
		}
	
		cv::Mat temp;

		//返回结果AL
		temp = GetAandL(lpt, pt3d, i);

		for (int j = 0; j < temp.rows; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				AL.at<float>(j + jsq, k) = temp.at<float>(j, k);
			}
		}
		jsq = jsq + temp.rows;
	}

	return AL;

}



//返回双/单像的AL矩阵，若是第一张和第二张的，返回两张的AL，除此之外都是返回右片的AL,flag表达是否第一张和第二张，
cv::Mat CMulti_SpaceIntersection::GetAandL(SATPoint2D* lpt, SATPoint3D* ptObj, int flag)
{
	//同名点数量
	int npt = 1;
	//flag是上面的循环次，将第一次和后面的次分开

	double A[2][3], l[2][1];

	double Pn = 0, Ln = 0, Hn = 0;
	double NumL = 0, DenL = 0, NumS = 0, DenS = 0;
	double dNumLdLn = 0.0, dNumLdPn = 0.0, dNumLdHn = 0.0;
	double dDenLdLn = 0.0, dDenLdPn = 0.0, dDenLdHn = 0.0;
	double dNumSdLn = 0.0, dNumSdPn = 0.0, dNumSdHn = 0.0;
	double dDenSdLn = 0.0, dDenSdPn = 0.0, dDenSdHn = 0.0;  //上面四行用来记录正解法偏导数

	cv::Mat AandL(2 * npt, 4, CV_32F);
	double Ltmp, Stmp;

	for (int i = 0; i < npt; i++)
	{
		// while(numofiterative<30&&(fabs(d[0][0])>0.00001||fabs(d[1][0])>0.00001||fabs(d[2][0])>0.001))
		//……………………处理左片………………	
		Pn = (ptObj[i].P - p_lRPC.LAT_OFF) / p_lRPC.LAT_SCALE;
		Ln = (ptObj[i].L - p_lRPC.LONG_OFF) / p_lRPC.LONG_SCALE;
		Hn = (ptObj[i].H - p_lRPC.HEIGHT_OFF) / p_lRPC.HEIGHT_SCALE;

		NumL = getaccumulation1(p_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
		DenL = getaccumulation1(p_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
		NumS = getaccumulation1(p_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
		DenS = getaccumulation1(p_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);

		dNumLdLn = getpartialderivativeofL1(p_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
		dNumLdPn = getpartialderivativeofP1(p_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
		dNumLdHn = getpartialderivativeofH1(p_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);

		dDenLdLn = getpartialderivativeofL1(p_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
		dDenLdPn = getpartialderivativeofP1(p_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
		dDenLdHn = getpartialderivativeofH1(p_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);

		dNumSdLn = getpartialderivativeofL1(p_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
		dNumSdPn = getpartialderivativeofP1(p_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
		dNumSdHn = getpartialderivativeofH1(p_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);

		dDenSdLn = getpartialderivativeofL1(p_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);
		dDenSdPn = getpartialderivativeofP1(p_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);
		dDenSdHn = getpartialderivativeofH1(p_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);

		A[0][0] = (p_lRPC.LINE_SCALE / p_lRPC.LONG_SCALE) * ((dNumLdLn * DenL - dDenLdLn * NumL) / (DenL * DenL));
		A[0][1] = (p_lRPC.LINE_SCALE / p_lRPC.LAT_SCALE) * ((dNumLdPn * DenL - dDenLdPn * NumL) / (DenL * DenL));
		A[0][2] = (p_lRPC.LINE_SCALE / p_lRPC.HEIGHT_SCALE) * ((dNumLdHn * DenL - dDenLdHn * NumL) / (DenL * DenL));
		A[1][0] = (p_lRPC.SAMP_SCALE / p_lRPC.LONG_SCALE) * ((dNumSdLn * DenS - dDenSdLn * NumS) / (DenS * DenS));
		A[1][1] = (p_lRPC.SAMP_SCALE / p_lRPC.LAT_SCALE) * ((dNumSdPn * DenS - dDenSdPn * NumS) / (DenS * DenS));
		A[1][2] = (p_lRPC.SAMP_SCALE / p_lRPC.HEIGHT_SCALE) * ((dNumSdHn * DenS - dDenSdHn * NumS) / (DenS * DenS));

		Ltmp = ((NumL / DenL) * p_lRPC.LINE_SCALE + p_lRPC.LINE_OFF);
		Stmp = ((NumS / DenS) * p_lRPC.SAMP_SCALE + p_lRPC.SAMP_OFF);

		double L = p_lAffine.lineb0 + p_lAffine.lineb1 * Stmp + p_lAffine.lineb2 * Ltmp;
		double S = p_lAffine.samplea0 + p_lAffine.samplea1 * Stmp + p_lAffine.samplea2 * Ltmp;

		l[0][0] = lpt[i].line - L;
		l[1][0] = lpt[i].sample - S;



		//……………………处理右片………………………………


		for (int i1 = 0; i1 < 2; i1++)
		{
			for (int j1 = 0; j1 < 3; j1++)
			{
				AandL.at<float>(i1, j1) = A[i1][j1];
			}
			AandL.at<float>(i1, 3) = l[i1][0];
		}

		return AandL;
	}

}

//前方交会
bool CMulti_SpaceIntersection::RPCIntersection(SATPoint2D* lpt, SATPoint2D* rpt, int npt, SATPoint3D* ptObj)
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

	for (int i = 0; i < npt; i++)
	{
		ptObj[i].L = (p_lRPC.LONG_OFF + p_rRPC.LONG_OFF) / 2;
		ptObj[i].P = (p_lRPC.LAT_OFF + p_rRPC.LAT_OFF) / 2;
		ptObj[i].H = (p_lRPC.HEIGHT_OFF + p_rRPC.HEIGHT_OFF) / 2;
	}

	for (int i = 0; i < npt; i++)
	{
		// while(numofiterative<30&&(fabs(d[0][0])>0.00001||fabs(d[1][0])>0.00001||fabs(d[2][0])>0.001))
		while (numofiterative < 15)
		{
			//处理左片

			Pn = (ptObj[i].P - p_lRPC.LAT_OFF) / p_lRPC.LAT_SCALE;
			Ln = (ptObj[i].L - p_lRPC.LONG_OFF) / p_lRPC.LONG_SCALE;
			Hn = (ptObj[i].H - p_lRPC.HEIGHT_OFF) / p_lRPC.HEIGHT_SCALE;

			NumL = getaccumulation1(p_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
			DenL = getaccumulation1(p_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
			NumS = getaccumulation1(p_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
			DenS = getaccumulation1(p_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);

			dNumLdLn = getpartialderivativeofL1(p_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
			dNumLdPn = getpartialderivativeofP1(p_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
			dNumLdHn = getpartialderivativeofH1(p_lRPC.LINE_NUM_COEFF, Ln, Pn, Hn);

			dDenLdLn = getpartialderivativeofL1(p_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
			dDenLdPn = getpartialderivativeofP1(p_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
			dDenLdHn = getpartialderivativeofH1(p_lRPC.LINE_DEN_COEFF, Ln, Pn, Hn);

			dNumSdLn = getpartialderivativeofL1(p_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
			dNumSdPn = getpartialderivativeofP1(p_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
			dNumSdHn = getpartialderivativeofH1(p_lRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);

			dDenSdLn = getpartialderivativeofL1(p_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);
			dDenSdPn = getpartialderivativeofP1(p_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);
			dDenSdHn = getpartialderivativeofH1(p_lRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);

			A[0][0] = (p_lRPC.LINE_SCALE / p_lRPC.LONG_SCALE) * ((dNumLdLn * DenL - dDenLdLn * NumL) / (DenL * DenL));
			A[0][1] = (p_lRPC.LINE_SCALE / p_lRPC.LAT_SCALE) * ((dNumLdPn * DenL - dDenLdPn * NumL) / (DenL * DenL));
			A[0][2] = (p_lRPC.LINE_SCALE / p_lRPC.HEIGHT_SCALE) * ((dNumLdHn * DenL - dDenLdHn * NumL) / (DenL * DenL));
			A[1][0] = (p_lRPC.SAMP_SCALE / p_lRPC.LONG_SCALE) * ((dNumSdLn * DenS - dDenSdLn * NumS) / (DenS * DenS));
			A[1][1] = (p_lRPC.SAMP_SCALE / p_lRPC.LAT_SCALE) * ((dNumSdPn * DenS - dDenSdPn * NumS) / (DenS * DenS));
			A[1][2] = (p_lRPC.SAMP_SCALE / p_lRPC.HEIGHT_SCALE) * ((dNumSdHn * DenS - dDenSdHn * NumS) / (DenS * DenS));

			Ltmp = ((NumL / DenL) * p_lRPC.LINE_SCALE + p_lRPC.LINE_OFF);
			Stmp = ((NumS / DenS) * p_lRPC.SAMP_SCALE + p_lRPC.SAMP_OFF);

			double Ll = p_lAffine.lineb0 + p_lAffine.lineb1 * Stmp + p_lAffine.lineb2 * Ltmp;
			double Sl = p_lAffine.samplea0 + p_lAffine.samplea1 * Stmp + p_lAffine.samplea2 * Ltmp;

			l[0][0] = lpt[i].line - Ll;
			l[1][0] = lpt[i].sample - Sl;

			//……………………处理右片………………………………

			Pn = (ptObj[i].P - p_rRPC.LAT_OFF) / p_rRPC.LAT_SCALE;
			Ln = (ptObj[i].L - p_rRPC.LONG_OFF) / p_rRPC.LONG_SCALE;
			Hn = (ptObj[i].H - p_rRPC.HEIGHT_OFF) / p_rRPC.HEIGHT_SCALE;


			NumL = getaccumulation1(p_rRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
			DenL = getaccumulation1(p_rRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
			NumS = getaccumulation1(p_rRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
			DenS = getaccumulation1(p_rRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);

			dNumLdLn = getpartialderivativeofL1(p_rRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
			dNumLdPn = getpartialderivativeofP1(p_rRPC.LINE_NUM_COEFF, Ln, Pn, Hn);
			dNumLdHn = getpartialderivativeofH1(p_rRPC.LINE_NUM_COEFF, Ln, Pn, Hn);

			dDenLdLn = getpartialderivativeofL1(p_rRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
			dDenLdPn = getpartialderivativeofP1(p_rRPC.LINE_DEN_COEFF, Ln, Pn, Hn);
			dDenLdHn = getpartialderivativeofH1(p_rRPC.LINE_DEN_COEFF, Ln, Pn, Hn);

			dNumSdLn = getpartialderivativeofL1(p_rRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
			dNumSdPn = getpartialderivativeofP1(p_rRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);
			dNumSdHn = getpartialderivativeofH1(p_rRPC.SAMP_NUM_COEFF, Ln, Pn, Hn);

			dDenSdLn = getpartialderivativeofL1(p_rRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);
			dDenSdPn = getpartialderivativeofP1(p_rRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);
			dDenSdHn = getpartialderivativeofH1(p_rRPC.SAMP_DEN_COEFF, Ln, Pn, Hn);

			A[2][0] = (p_rRPC.LINE_SCALE / p_rRPC.LONG_SCALE) * ((dNumLdLn * DenL - dDenLdLn * NumL) / (DenL * DenL));
			A[2][1] = (p_rRPC.LINE_SCALE / p_rRPC.LAT_SCALE) * ((dNumLdPn * DenL - dDenLdPn * NumL) / (DenL * DenL));
			A[2][2] = (p_rRPC.LINE_SCALE / p_rRPC.HEIGHT_SCALE) * ((dNumLdHn * DenL - dDenLdHn * NumL) / (DenL * DenL));
			A[3][0] = (p_rRPC.SAMP_SCALE / p_rRPC.LONG_SCALE) * ((dNumSdLn * DenS - dDenSdLn * NumS) / (DenS * DenS));
			A[3][1] = (p_rRPC.SAMP_SCALE / p_rRPC.LAT_SCALE) * ((dNumSdPn * DenS - dDenSdPn * NumS) / (DenS * DenS));
			A[3][2] = (p_rRPC.SAMP_SCALE / p_rRPC.HEIGHT_SCALE) * ((dNumSdHn * DenS - dDenSdHn * NumS) / (DenS * DenS));

			Ltmp = ((NumL / DenL) * p_rRPC.LINE_SCALE + p_rRPC.LINE_OFF);
			Stmp = ((NumS / DenS) * p_rRPC.SAMP_SCALE + p_rRPC.SAMP_OFF);

			double Lr = p_rAffine.lineb0 + p_rAffine.lineb1 * Stmp + p_rAffine.lineb2 * Ltmp;
			double Sr = p_rAffine.samplea0 + p_rAffine.samplea1 * Stmp + p_rAffine.samplea2 * Ltmp;

			l[2][0] = rpt[i].line - Lr;
			l[3][0] = rpt[i].sample - Sr;

			for (int p = 0; p < 4; p++)
			{
				for (int q = 0; q < 3; q++)
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
			else
			{
				// std::cout << rpt[i].line <<"  " << rpt[i].sample <<"  "<< lpt[i].line << "  "<< lpt[i].sample << std::endl;
				// std::cout << Lr << "  " << Sr << "  " << Ltmp << "  " << Stmp << std::endl;
				// cout<<endl<<"timestimes\t"<<timestimes<<endl;
				return false;
			}

			ptObj[i].L += d[0][0];

			ptObj[i].P += d[1][0];

			ptObj[i].H += d[2][0];

			numofiterative++;

		}
		numofiterative = 0;
		d[0][0] = 1;
		d[1][0] = 1;
		d[2][0] = 50;

	}
	return true;
	
}

//获取仿射变换的逆矩阵
// 	L0=lineb0+lineb1*L+lineb2*S      //右边地面坐标的投影
//  S0=samplea0+samplea1*S+samplea2*L //原始影像上的坐标
//  L = dstlineb0+dstlineb1*L0+dstlineb2*S0
//  S = dstsamplea0+dstsamplea1*L0+dstsamplea2*s0
bool CMulti_SpaceIntersection::GetInverseAffPara(RPCImAffine& srcAffPara, RPCImAffine& dstAffPara)
{
	double b0 = srcAffPara.lineb0, b1 = srcAffPara.lineb1, b2 = srcAffPara.lineb2;
	double a0 = srcAffPara.samplea0, a1 = srcAffPara.samplea1, a2 = srcAffPara.samplea2;

	double fenmu = a1 * b2 - a2 * b1;

	if (fenmu == 0)
	{
		return false;
	}
	//
	//
	double inverA[6];
	//Line 相关
	inverA[0] = (a2 * b0 - a0 * b2) / fenmu;
	inverA[1] = b2 / fenmu;
	inverA[2] = a2 / fenmu;


	//SAMPLE 相关
	inverA[3] = (a0 * b1 - a1 * b0) / fenmu;
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

//给定一个高程，以及对应的影像RPC参数，计算影像（sample，line）对应的 地面坐标
void CMulti_SpaceIntersection::RPCImg2Obj(RPCcoeffcient& RPCcoef, double H, RPCImAffine& affpara, SATPoint2D pimgpt, SATPoint3D& ObjPt)
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

//空间三位点投影到影像
SATPoint2D CMulti_SpaceIntersection::RPCObj2Img(RPCcoeffcient& RPCcoef, SATPoint3D& ObjPt, RPCImAffine& affpara)
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
	return ptImg;
}

//求逆矩阵
int CMulti_SpaceIntersection::inv(double* m1, int n)
{
	int* is, * js;
	int i, j, k, l, u, v;
	double temp, max_v;
	is = (int*)malloc(n * sizeof(int));
	js = (int*)malloc(n * sizeof(int));
	if (is == NULL || js == NULL)
	{
		printf("out of memory!\n");
		return(0);
	}
	for (k = 0; k < n; k++)
	{
		max_v = 0.0;
		for (i = k; i < n; i++)
			for (j = k; j < n; j++)
			{
				temp = fabs(m1[i * n + j]);
				if (temp > max_v)
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
			for (j = 0; j < n; j++)
			{
				u = k * n + j; v = is[k] * n + j;
				temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
			}
		if (js[k] != k)
			for (i = 0; i < n; i++)
			{
				u = i * n + k; v = i * n + js[k];
				temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
			}
		l = k * n + k;
		m1[l] = 1.0 / m1[l];
		for (j = 0; j < n; j++)
			if (j != k)
			{
				u = k * n + j;
				m1[u] *= m1[l];
			}
		for (i = 0; i < n; i++)
			if (i != k)
				for (j = 0; j < n; j++)
					if (j != k)
					{
						u = i * n + j;
						m1[u] -= m1[i * n + k] * m1[k * n + j];
					}
		for (i = 0; i < n; i++)
			if (i != k)
			{
				u = i * n + k;
				m1[u] *= -m1[l];
			}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j < n; j++)
			{
				u = k * n + j; v = js[k] * n + j;
				temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
			}
		if (is[k] != k)
			for (i = 0; i < n; i++)
			{
				u = i * n + k; v = i * n + is[k];
				temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
			}
	}
	free(is);
	free(js);
	return(1);
}

//矩阵相乘
void CMulti_SpaceIntersection::matrixmulti(double* r1, double* r2, double* r, int m, int n, int p)
{
	for (int i = 0; i < m * p; i++)
		*(r + i) = 0;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < p; j++)
		{
			for (int k = 0; k < n; k++)
			{
				*(r + i * p + j) += *(r1 + i * n + k) * (*(r2 + k * p + j));
			}

		}
	}
}



//求偏导
double CMulti_SpaceIntersection::getpartialderivativeofL1(double* rpc, double L, double P, double H)
{
	/*cout << setprecision(44) << rpc[1] << endl;*/
	return(rpc[1] + rpc[4] * P + rpc[5] * H + 2 * rpc[7] * L + rpc[10] * P * H + 3 * rpc[11] * L * L
		+ rpc[12] * P * P + rpc[13] * H * H + 2 * rpc[14] * P * L + 2 * rpc[17] * L * H);
}
double CMulti_SpaceIntersection::getpartialderivativeofP1(double* rpc, double L, double P, double H)
{
	return(rpc[2] + rpc[4] * L + rpc[6] * H + 2 * rpc[8] * P + rpc[10] * L * H + 2 * rpc[12] * L * P
		+ rpc[14] * L * L + 3 * rpc[15] * P * P + rpc[16] * H * H + 2 * rpc[18] * P * H);
}
double CMulti_SpaceIntersection::getpartialderivativeofH1(double* rpc, double L, double P, double H)
{
	return(rpc[3] + rpc[5] * L + rpc[6] * P + 2 * rpc[9] * H + rpc[10] * P * L + 2 * rpc[13] * L * H
		+ 2 * rpc[16] * P * H + rpc[17] * L * L + rpc[18] * P * P + 3 * rpc[19] * H * H);
}
//求和
double CMulti_SpaceIntersection::getaccumulation1(double* rpc, double L, double P, double H)
{
	double M[20];
	double S = 0;
	M[0] = 1, M[1] = L, M[2] = P, M[3] = H,
		M[4] = L * P, M[5] = L * H, M[6] = P * H, M[7] = L * L, M[8] = P * P, M[9] = H * H,
		M[10] = P * L * H, M[11] = L * L * L, M[12] = L * P * P, M[13] = L * H * H, M[14] = L * L * P,
		M[15] = P * P * P, M[16] = P * H * H, M[17] = L * L * H, M[18] = P * P * H, M[19] = H * H * H;

	for (int j = 0; j < 20; j++)
	{
		S += rpc[j] * M[j];
	}
	return S;

}
//求和
double CMulti_SpaceIntersection::getaccumulation1(double rpc[20], SATPoint3D& ObjPt)
{
	double M[20];
	double S = 0;
	double L = ObjPt.L;
	double P = ObjPt.P;
	double H = ObjPt.H;
	M[0] = 1, M[1] = L, M[2] = P, M[3] = H,
		M[4] = L * P, M[5] = L * H, M[6] = P * H, M[7] = L * L, M[8] = P * P, M[9] = H * H,
		M[10] = P * L * H, M[11] = L * L * L, M[12] = L * P * P, M[13] = L * H * H, M[14] = L * L * P,
		M[15] = P * P * P, M[16] = P * H * H, M[17] = L * L * H, M[18] = P * P * H, M[19] = H * H * H;

	for (int j = 0; j < 20; j++)
	{
		S += rpc[j] * M[j];
	}
	return S;

}
//求与P相关的偏导数
double CMulti_SpaceIntersection::GetPartialDerivativeofP(double Numrpc[20], double Denrpc[20], SATPoint3D& objpt, double SL)
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
double CMulti_SpaceIntersection::GetPartialDerivativeofL(double Numrpc[20], double Denrpc[20], SATPoint3D& objpt, double SL)
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