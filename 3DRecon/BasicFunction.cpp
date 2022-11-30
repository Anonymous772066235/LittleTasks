#include "BasicFunction.h"
#include <stdio.h>

using namespace std;

//去除文件中空格
void BasicFunction::Trim(std::string & str)
{
	std::string blanks("\f\v\r\t\n ");
    str.erase(0,str.find_first_not_of(blanks));
    str.erase(str.find_last_not_of(blanks) + 1);
}
//读取上级文件夹目录 
bool BasicFunction::GetFolderPath(std::string FilePath, std::string &FolderPath)
{
	auto Iter = std::find(FilePath.crbegin(), FilePath.crend(), '/');
	FolderPath = std::string(FilePath.cbegin(), Iter.base());
	return 1;
	
}
//按行读取txt中路径
bool BasicFunction::ReadImageNames(std::string ListPath, vector<std::string> &Data)
{
	ifstream infile;

	infile.open(ListPath);   //将文件流对象与文件连接起来 
	if (infile.is_open())//若失败,则输出错误消息,并终止程序运行 
	{
		string s;
		while (getline(infile, s))
		{
			Trim(s);
			Data.push_back(s.c_str());
			std::cout<<s<<endl;
		}
	}
	else
	{
		return 0;
	}
	infile.close();             //关闭文件输入
	return 1;
}

//日志-----------------------------------------------------------------------------------------------
//创建文件夹
int  BasicFunction::CreateFolder(std::string FolderPath)
{
	int State = mkdir(FolderPath.c_str(),0);
	return State;
}

//删除文件夹
int BasicFunction::DeleteFolder(std::string FolderPath)
{
	int State = rmdir(FolderPath.c_str());
	return State;
}

//写入指定string到文件(没有则创建，不覆盖)
void BasicFunction::Write2Text(std::string FilePath, std::string Data)
{
	ofstream os;     //创建一个文件输出流对象
	os.open(FilePath, ios::out);//将对象与文件关联，ios::app 追加，不覆盖
	os << Data << endl;
	os.close();
}

//删除指定文件
int BasicFunction::DeleteFile(std::string FilePath)
{
	if (remove(FilePath.c_str()) == 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
//查询文件是否存在
bool BasicFunction::File_Exist(const std::string& FilePath)
{
	ifstream f(FilePath.c_str());
	return f.good();
}
bool BasicFunction::Path_Exist(const std::string& name)
{
	return access(name.c_str(), 0);
}
//将指定命令转化为json格式，并存到指定文件
std::string BasicFunction::Convert2Json_new(int Code, std::string Message, float Process, std::string isEnd)
{
	std::string JsonMessage;
	
	{
		//获取系统时间
		time_t timep;
		struct tm* p;
		time(&timep); //获取从1970至今过了多少秒，存入time_t类型的timep
		p = localtime(&timep);//用localtime将秒数转化为struct tm结构体

		//将时间转string
		std::string sTime = std::to_string(p->tm_year + 1900) + "-" + std::to_string(p->tm_mon+1) + "-" + std::to_string(p->tm_mday) +
			" " + std::to_string(p->tm_hour) + ":" + std::to_string(p->tm_min) + ":" + std::to_string(p->tm_sec);
		//数据
		//进度
		JsonMessage = "{\"process\":\"";
		JsonMessage += std::to_string(Process) + "\",";

		//时间
		JsonMessage += "\"time\":\"";
		JsonMessage += sTime;
		JsonMessage += "\",";


		//错误代码
		JsonMessage += "\"code\":\"";
		JsonMessage += std::to_string(Code) + "\",";

		//消息
		JsonMessage += "\"msg\":\"";
		JsonMessage += Message + "\",";

		//消息
		JsonMessage += "\"isEnd\":\"";
		JsonMessage += isEnd + "\"}";

		return JsonMessage;
	}
}

// original version of Convert2Json
string BasicFunction::Convert2Json(int Code, std::string Message, float Process, std::string CurrentWork)
{
	std::string JsonMessage;
	if (Process == 0)
	{
		//数据为空
		JsonMessage = "{\"data\":\"\",";

		//错误代码
		JsonMessage += "\"code\":";
		JsonMessage += std::to_string(Code) + ",";

		//消息
		JsonMessage += "\"msg\":\"";
		JsonMessage += Message + "\"}";

		return JsonMessage;
	}
	else
	{
		//获取系统时间
		time_t timep;
		struct tm* p;
		time(&timep); //获取从1970至今过了多少秒，存入time_t类型的timep
		p = localtime(&timep);//用localtime将秒数转化为struct tm结构体

		//将时间转string
		std::string sTime = std::to_string(p->tm_year + 1900) + "-" + std::to_string(p->tm_mon+1) + "-" + std::to_string(p->tm_mday) +
			" " + std::to_string(p->tm_hour) + ":" + std::to_string(p->tm_min) + ":" + std::to_string(p->tm_sec);
		//数据
		//进度
		JsonMessage = "{\"data\":[{\"process\":";
		JsonMessage += std::to_string(Process) + ",";

		//当前工作
		JsonMessage += "\"name\":\"";
		JsonMessage += CurrentWork + "\",";

		//时间
		JsonMessage += "\"time\":\"";
		JsonMessage += sTime;
		JsonMessage += "\"}],";


		//错误代码
		JsonMessage += "\"code\":";
		JsonMessage += std::to_string(Code) + ",";

		//消息
		JsonMessage += "\"msg\":\"";
		JsonMessage += Message + "\"}";

		return JsonMessage;
	}
}
//读取指定后缀名的文件
void BasicFunction::LoadFilesPath(string FolderPath, string Exd, vector<string>& FilesPath)
{
	/*const std::string path0 = FolderPath;
    DIR* pDir;
    struct dirent* ptr;
 
    struct stat s;
    lstat(FolderPath.c_str(), &s);
 
    if(!S_ISDIR(s.st_mode)){
        std::cout << "not a valid directory: " << FolderPath << std::endl;
        return;
    }
 
    if(!(pDir = opendir(FolderPath.c_str()))){
        std::cout << "opendir error: " << FolderPath << std::endl;
        return;
    }
    int i = 0;
    std::string subFile;
    while((ptr = readdir(pDir)) != 0){
        subFile = ptr -> d_name;
        if(subFile == "." || subFile == "..")
            continue;
		if(BasicFunction::GetFileSuffix(subFile)==Exd)
		{
        subFile = path0 + subFile;
        std::cout << ++i << ": " << subFile << std::endl;
        FilesPath.push_back(subFile);
		}
    }
    closedir(pDir);*/

}


cv::Mat BasicFunction::ReadTifoneBand(std::string filename)
{
	char* Proj;
	double* GeoTransform;
	//初始化GDAL库注册表
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	//读取影像
	GDALDataset* ImageData = (GDALDataset*)GDALOpen(filename.data(), GA_ReadOnly);
	//获取投影转换
	const char* const_prj = ImageData->GetProjectionRef();
	Proj = const_cast<char*>(const_prj);
	//std::cout << Proj << std::endl;
	//获取地理信息
	GeoTransform = new double[6];
	ImageData->GetGeoTransform(GeoTransform);
	/*for (int i(0); i<6; ++i)
	{
		cout << GeoTransform[i] << endl;
	}*/
	//获取栅格数据，默认只保留第一波段
	int Width = ImageData->GetRasterXSize(); // 列 
	int Height = ImageData->GetRasterYSize(); // 行 
	//int BandSize = ImageData->GetRasterCount();//波段数
	float* Data_of_Band;   // 存储数据
	Data_of_Band = new float[Width * Height];
	GDALRasterBand* pBand = ImageData->GetRasterBand(1);
	GDALDataType DataType = pBand->GetRasterDataType();
	pBand->RasterIO(GF_Read, 0, 0, Width, Height, Data_of_Band,
		Width, Height, GDT_Float32, 0, 0);
	cv::Mat Band_Image = cv::Mat(Height, Width, CV_32FC1, Data_of_Band);
	//2%线性拉伸的方式转换为uint8
	cv::Mat nimg = LinearStrech(Band_Image, 0.01);
	Band_Image.release();
	if (Data_of_Band != NULL)
	{
		delete[] Data_of_Band;
	}
	return nimg;
}


std::vector<cv::Mat> BasicFunction::ReadTif(std::string filename)
{
	char* Proj;
	double* GeoTransform;
	//初始化GDAL库注册表
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	//读取影像
	GDALDataset* ImageData = (GDALDataset*)GDALOpen(filename.data(), GA_ReadOnly);
	//获取投影转换
	const char* const_prj = ImageData->GetProjectionRef();
	Proj = const_cast<char*>(const_prj);
	//获取地理信息
	GeoTransform = new double[6];
	ImageData->GetGeoTransform(GeoTransform);
	//获取栅格数据，默认只保留第一波段
	int Width = ImageData->GetRasterXSize(); // 列 
	int Height = ImageData->GetRasterYSize(); // 行 
	int BandSize = ImageData->GetRasterCount();//波段数
	std::vector<cv::Mat> imgs;
	for (int i(0); i < BandSize; ++i)
	{
		float* Data_of_Band;   // 存储数据
		Data_of_Band = new float[Width * Height];
		GDALRasterBand* pBand = ImageData->GetRasterBand(i + 1);
		GDALDataType DataType = pBand->GetRasterDataType();
		pBand->RasterIO(GF_Read, 0, 0, Width, Height, Data_of_Band,
			Width, Height, GDT_Float32, 0, 0);
		cv::Mat Band_Image = cv::Mat(Height, Width, CV_32FC1, Data_of_Band);
		//2%线性拉伸的方式转换为uint8
		cv::Mat nimg = LinearStrech(Band_Image, 0.01);
		Band_Image.release();
		if (Data_of_Band != NULL)
		{
			delete[] Data_of_Band;
		}
		imgs.push_back(nimg);
	}
	return imgs;
}


bool BasicFunction::GetBoundaryfromXML(std::string xmlpath, Boundary& bound)
{
	tinyxml2::XMLDocument XmlFile;
	if (XmlFile.LoadFile(xmlpath.c_str()) != tinyxml2::XML_SUCCESS)
	{
		std::cout << "打开xml文件失败，请检查文件目录!" << std::endl;
		return false;
	}
	//tinyxml2::XMLElement* rootElement = XmlFile.RootElement()->FirstChildElement("ProductMetaData");
	tinyxml2::XMLElement* rootElement = XmlFile.RootElement();
	const char* tMetadata;
	//左上角
	if (rootElement->FirstChildElement("TopLeftLongitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("TopLeftLongitude")->GetText();
		bound.TopLeftLongitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取左上角经度失败";
	}
	if (rootElement->FirstChildElement("TopLeftLatitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("TopLeftLatitude")->GetText();
		bound.TopLeftLatitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取左上角纬度度失败";
	}
	//右上角
	if (rootElement->FirstChildElement("TopRightLongitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("TopRightLongitude")->GetText();
		bound.TopRightLongitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取右上角经度失败";
	}
	if (rootElement->FirstChildElement("TopRightLatitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("TopRightLatitude")->GetText();
		bound.TopRightLatitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取右上角经度失败";
	}
	//左下角
	if (rootElement->FirstChildElement("BottomLeftLongitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("BottomLeftLongitude")->GetText();
		bound.BottomLeftLongitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取左下角经度失败";
	}
	if (rootElement->FirstChildElement("BottomLeftLatitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("BottomLeftLatitude")->GetText();
		bound.BottomLeftLatitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取左下角经度失败";
	}
	//右下角
	if (rootElement->FirstChildElement("BottomRightLongitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("BottomRightLongitude")->GetText();
		bound.BottomRightLongitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取左上角经度失败";
	}
	if (rootElement->FirstChildElement("BottomRightLatitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("BottomRightLatitude")->GetText();
		bound.BottomRightLatitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取左上角经度失败";
	}
	return true;
}

//读取多波段影像（多个文件）
ImageAndGeoinfo BasicFunction::ReadMulBandImageFromFolder(std::string ImagePath, GDALRPCInfo &RPCParams)
{
	//首先读取指定影像
	GDALDataset* FirstBand = BasicFunction::ReadGTiffImage(ImagePath, RPCParams);

	//指定影像的尺寸
	int Height = FirstBand->GetRasterYSize();
	int Width = FirstBand->GetRasterXSize();

	//分割文件名和文件夹路径
	std::string FolderPath, FileName;
	BasicFunction::SegFilePath(ImagePath, FolderPath, FileName);

	//获取当前文件的后缀名
	std::string Suffix = BasicFunction::GetFileSuffix(ImagePath);
	std::vector<std::string>ImageFilePath;
	BasicFunction::LoadFilesPath(FolderPath, Suffix, ImageFilePath);
	

	//需要返回的额信息
	ImageAndGeoinfo ImageData;
	//读取信息
	for (int i = 0; i < ImageFilePath.size(); i++)
	{
		GDALRPCInfo iRPCParams;
		GDALDataset* iImageData = BasicFunction::ReadGTiffImage(ImageFilePath[i], iRPCParams);
		if (iImageData->GetRasterYSize() == Height && iImageData->GetRasterXSize() == Width)
		{
			//读取每个波段
			std::vector<cv::Mat> iBand;

			//
			char* iProjParams = new char[1000];
			double iGeoParams[6];
			if (BasicFunction::GDAL2Mat(iImageData, iBand, iProjParams, iGeoParams))
			{
				for (int j = 0; j < iBand.size(); j++)
				{
					ImageData._Image.push_back(iBand[j]);
				}
			}
			if (i == 0)
			{
				ImageData._GeoTrans = iGeoParams;
				ImageData._Proj = iProjParams;
			}
		}
		//iImageData->Release();
	}

	return ImageData;
}
//读取GF2影像及RPB文件
GDALDataset* BasicFunction::ReadGTiffImage(std::string ImagePath, GDALRPCInfo &RPCParams)
{
	//初始化GDAL库注册表
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	//分割文件名和文件夹路径
	std::string FolderPath, FileName;
	BasicFunction::SegFilePath(ImagePath, FolderPath, FileName);

	//获取当前文件的后缀名
	std::vector<std::string>ImageFilePath;
	BasicFunction::LoadFilesPath(FolderPath, "", ImageFilePath);

	for (int i = 0; i < ImageFilePath.size(); i++)
	{
		string iSuffix = BasicFunction::GetFileSuffix(ImageFilePath[i]);
		if (iSuffix == "rpb")
		{
			break;
		}
		else if (iSuffix == "rpc")
		{
			string NewImageFilePath;
			ChangeSuffix(ImageFilePath[i], "rpb", NewImageFilePath);
			int ret = rename(ImageFilePath[i].data(), NewImageFilePath.data());

			//break;
		}

	}

	//读取影像
	GDALDataset* ImageData = (GDALDataset*)GDALOpen(ImagePath.data(), GA_ReadOnly);
	
	if (ImageData == NULL)
	{
		std::cout << "读取影像数据失败，请检查影像目录" << std::endl;
		return FALSE;
	}
	else
	{

		//读取rpc参数
		char** RPC = ImageData->GetMetadata("RPC");
		GDALExtractRPCInfo(RPC, &RPCParams);//格式化到结构体oInfo变量中
		
		delete[]RPC;
		return ImageData;
	}
}
//读取hdf影像
GDALDataset* BasicFunction::ReadHDFImage(std::string ImagePath)
{

	//初始化GDAL库注册表
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	//读取影像
	GDALDataset* ImageData_Total = (GDALDataset*)GDALOpen(ImagePath.data(), GA_ReadOnly);

	vector<string> vtDatasets;
	char** metadata = ImageData_Total->GetMetadata("SUBDATASETS");
	if (0 < CSLCount(const_cast<const char* const*>(metadata)))
	{
		for (int i = 0; metadata[i]; i ++)
		{
			string strDSName(metadata[i]);
			//std::cout << strDSName << '\n';
			strDSName = strDSName.substr(strDSName.find_first_of('=') + 1);
			vtDatasets.push_back(strDSName);
		}
	}

	GDALDataset* ImageData =  (GDALDataset*)GDALOpen(vtDatasets[4].c_str(), GA_ReadOnly);
	return ImageData;

}

bool isPointInRect(cv::Point p, cv::Rect rect)
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
//
// 鲁棒的surfcuda特征匹配
void BasicFunction::RobotSURF_CUDAMatch(const cv::Mat img1, const cv::Mat img2,
	vector<cv::Point2f>& p1, vector<cv::Point2f>& p2)
{
	//特征提取
	cv::cuda::GpuMat keypoints1, descriptions1, keypoints2, descriptions2;
	cv::cuda::GpuMat dimg1(img1), dimg2(img2);
	cv::cuda::SURF_CUDA surf_cuda;
	surf_cuda(dimg1, cv::cuda::GpuMat(), keypoints1, descriptions1);
	surf_cuda(dimg2, cv::cuda::GpuMat(), keypoints2, descriptions2);
	//download关键点
	vector<cv::KeyPoint> kp1, kp2;
	surf_cuda.downloadKeypoints(keypoints1, kp1);
	surf_cuda.downloadKeypoints(keypoints2, kp2);
	//硬匹配
	auto matcher = cv::cuda::DescriptorMatcher::createBFMatcher(cv::NORM_L2);
	vector<cv::DMatch> matches;
	matcher->match(descriptions1, descriptions2, matches);
	//GMS剔除误匹配
	vector<cv::DMatch> mts;
	cv::xfeatures2d::matchGMS(img1.size(), img2.size(), kp1, kp2, matches, mts, 0, 0, 10.);
	//ransac剔除误匹配
	vector<cv::DMatch> GoodMatches = mts;
	//vector<cv::KeyPoint> RAN_KP1, RAN_KP2;
	//for (int i(0); i < GoodMatches.size(); ++i)
	//{
	//	RAN_KP1.push_back(kp1[mts[i].queryIdx]);
	//	RAN_KP2.push_back(kp2[mts[i].trainIdx]);
	//}
	////坐标数据类型变换
	//vector<cv::Point2f> p01, p02;
	//for (int j(0); j < GoodMatches.size(); ++j)
	//{
	//	p01.push_back(RAN_KP1[j].pt);
	//	p02.push_back(RAN_KP2[j].pt);
	//}
	////RANSAC过程
	//vector<uchar> RansacStatus;
	//cv::Mat Fundamental = cv::findFundamentalMat(p01, p02, RansacStatus, cv::FM_RANSAC);
	//int iter = 0;
	//for (int i(0); i < mts.size(); ++i)
	//{
	//	if (RansacStatus[i] == 0)
	//		GoodMatches.erase(GoodMatches.begin() + iter);
	//	else
	//		iter++;
	//}
	//整理点集
	for (int i(0); i < GoodMatches.size(); ++i)
	{
		int index1 = GoodMatches[i].queryIdx;
		int index2 = GoodMatches[i].trainIdx;
		p1.push_back(kp1[index1].pt);
		p2.push_back(kp2[index2].pt);
	}
}


void BasicFunction::RobotSURF_CUDAMatch(const cv::Mat img1, const cv::Mat img2, std::vector<cv::Point2f>& p1,
	std::vector<cv::Point2f>& p2, int patchheight, int patchwidth)
{
	//分块进行关键点提取
	int rows = img1.rows;
	int cols = img1.cols;
	//获取块数目
	int rpatchnum = (float)rows / patchheight;
	int cpatchnum = (float)cols / patchwidth;
	//重算块大小
	patchheight = (float)rows / rpatchnum;
	patchwidth = (float)cols / cpatchnum;
	for (int i(0); i < rpatchnum; ++i)
	{
		for (int j(0); j < cpatchnum; ++j)
		{
			//获取块内影像
			cv::Mat mimg1 = img1(cv::Rect(j * patchwidth, i * patchheight, patchwidth, patchheight)).clone();
			cv::Mat mimg2 = img2(cv::Rect(j * patchwidth, i * patchheight, patchwidth, patchheight)).clone();
			//块内影像匹配
			std::vector<cv::Point2f> mp1, mp2;
			RobotSURF_CUDAMatch(mimg1, mimg2, mp1, mp2);
			//恢复关键点坐标并入总栈
			for (int k(0); k < mp1.size(); ++k)
			{
				mp1[k].x += j * patchwidth;
				mp1[k].y += i * patchheight;
				mp2[k].x += j * patchwidth;
				mp2[k].y += i * patchheight;
				p1.push_back(mp1[k]);
				p2.push_back(mp2[k]);
			}
		}
	}
}

//Mat转GDAL格式
GDALDataset* BasicFunction::Mat2GDAL(std::vector<cv::Mat>& IntputImage, char* Proj, double* GeoTransform)
{
	GDALDataset* tmp;
}
//获取dem的边界
DEMBoundary BasicFunction::GetDEMBoundary(Boundary pImageBoundary)
{

	DEMBoundary pDEMBoundary;
	//将影像边界写入cv矩阵
	cv::Mat Lon(4, 1, CV_32F);
	cv::Mat Lat(4, 1, CV_32F);
	Lon.at<float>(0, 0) = pImageBoundary.BottomLeftLongitude;
	Lon.at<float>(0, 1) = pImageBoundary.BottomRightLongitude;
	Lon.at<float>(0, 2) = pImageBoundary.TopLeftLongitude;
	Lon.at<float>(0, 3) = pImageBoundary.TopRightLongitude;

	Lat.at<float>(0, 0) = pImageBoundary.BottomLeftLatitude;
	Lat.at<float>(0, 1) = pImageBoundary.BottomRightLatitude;
	Lat.at<float>(0, 2) = pImageBoundary.TopLeftLatitude;
	Lat.at<float>(0, 3) = pImageBoundary.TopRightLatitude;

	//提取极值
	double MinLon, MinLat, MaxLon, MaxLat;

	cv::minMaxIdx(Lon, &MinLon, &MaxLon);
	cv::minMaxIdx(Lat, &MinLat, &MaxLat);
	if (MinLon < 0)
	{
		MinLon = 0;
	}
	if (MinLat < 0)
	{
		MinLat = 0;
	}

	pDEMBoundary.MinP = MinLat;
	pDEMBoundary.MinL = MinLon;

	pDEMBoundary.MaxP = MaxLat;
	pDEMBoundary.MaxL = MaxLon;

	//std::cout << "DEM边界: " << std::endl;
	//std::cout << "MinLon: " << pDEMBoundary.MinL << "  MinLat" << pDEMBoundary.MinP << std::endl;
	//std::cout << "MaxLon: " << pDEMBoundary.MaxL << "  MaxLat" << pDEMBoundary.MaxP << std::endl;
	return pDEMBoundary;
}
//获取文件后缀名
std::string BasicFunction::GetFileSuffix(std::string FilePath)
{
	std::string Suffix;
	for (int i = FilePath.size() - 1; i > 0; i--)
	{
		if (FilePath[i] == '.')
		{
			Suffix = FilePath.substr(i + 1);
			return Suffix;
		}
	}
	Suffix = FilePath;
	return Suffix;
}
//分开文件夹路径和后缀名
bool BasicFunction::SegFilePath(std::string FilePath, std::string& FolderPath, std::string& FileName)
{
	int index = FilePath.find_last_of("\\");
	if (index == -1)
	{
		index = FilePath.find_last_of("/");
	}
	//std::cout << index << std::endl;
	//文件夹
	FolderPath = FilePath.substr(0, index);
	//文件名
	FileName = FilePath.substr(index + 1, -1);


}
//更改文件后缀名
bool BasicFunction::ChangeSuffix(std::string FilePath, std::string Suffix,std::string & NFilePath)
{

	auto Iter = std::find(FilePath.crbegin(), FilePath.crend(), '.');
	NFilePath = std::string(FilePath.cbegin(), Iter.base());
	NFilePath += Suffix;

	return true;
}

//将GDAL数据写出为TIFF

void BasicFunction::SaveImageAsTiff(std::vector<cv::Mat> Image, std::string SaveFilePath, char* Proj, double* GeoTransform)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	const int BandSize =Image.size();
	const int Height = Image[0].rows;
	const int Width = Image[0].cols;

	
	GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	GDALDataset* ImageData;
	
	ImageData = poDriver->Create(SaveFilePath.data(),Width, Height, BandSize, GDT_Float32, NULL);

	const char* Projection = Proj;
	ImageData->SetGeoTransform(GeoTransform);
	ImageData->SetProjection(Projection);

	int PatchSize = 1000;


	if (Width > 10 * PatchSize && Height > 10 * PatchSize)
	{
		int WLops = Width / PatchSize + 1;
		int HLops = Height / PatchSize + 1;

		//cout << float(Width) / float(PatchSize) <<"  " << Width / PatchSize << endl;
		//闂備浇娉曢崰鎾剁懅婵犻潧鍊搁埀顒傚枑瑜把囨煙妞嬪骸鍘撮柡浣规崌瀵?剟濡堕崱妤婁紦闂備浇娉曢崰鎰板几婵犳艾绠?梺鍨?儐椤ρ囨⒑鐠恒劌鏋戦悷娆欑畵閺佸秴鐣濋崘鎯ф疂闂佸綊顥撻崗姗�寮?幘璇叉?闁靛牆妫楅?锟?1
		if (float(Width) / float(PatchSize) == Width / PatchSize)
		{
			WLops = WLops - 1;
		}
		if (float(Height) / float(PatchSize) ==  Height / PatchSize)
		{
			HLops = HLops - 1;
		}
		//cout << WLops << "  " << HLops << endl;
		for (int i = 0; i < BandSize; i++)
		{
			GDALRasterBand* pBand = ImageData->GetRasterBand(i + 1);

			for (int j = 0; j < HLops; j++)
			{
				for (int k = 0; k < WLops; k++)
				{
					if (j != HLops - 1 && k != WLops - 1)
					{
						float* Data_of_Band = new float[PatchSize * PatchSize];

						cv::Mat _kPatch;
						Image[i](cv::Rect(k* PatchSize, j* PatchSize, PatchSize, PatchSize)).copyTo(_kPatch );
						
						cv::Mat kPatch = _kPatch.reshape(1, 1);

						for (int l = 0; l < PatchSize * PatchSize; l++)
						{
							Data_of_Band[l] = kPatch.at<float>(0, l);

						}
						pBand->RasterIO(GF_Write, k * PatchSize, j * PatchSize, PatchSize, PatchSize, Data_of_Band,
							PatchSize, PatchSize, GDT_Float32, 0, 0);//闂備浇娉曢崰鎰板几婵犳艾绠?柣鎴ｅГ閺呮悂鏌￠崒妯衡偓鏍?偓姘鳖劑+1闂備浇娉曢崰鎰板几婵犳艾绠?柣鎴ｅГ閺呮悂鏌￠崒妯衡偓鏍?偓姘?秺閺屻劑鎮㈢拠娴嬫瀼缂佺偓婢橀ˇ鎵?偓姘?秺閺屻劑鎮㈤崨濠勪紕闂佸綊顥撻崗姗�寮?幘璇茬?闁汇値鍨辫ぐ褔鏌熸?搴″幋闁轰焦鎹囧?顒勫Χ閸℃?浼揇ata_of_Band

						delete[]Data_of_Band;
						_kPatch.release();
						kPatch.release();
					}
					else if (j != HLops - 1 && k == WLops - 1)
					{
						float* Data_of_Band = new float[PatchSize * (Width - k * PatchSize)];

						cv::Mat _kPatch;
						Image[i](cv::Rect(k* PatchSize, j* PatchSize, (Width - k * PatchSize), PatchSize)).copyTo(_kPatch);

						cv::Mat kPatch = _kPatch.reshape(1, 1);
						for (int l = 0; l < PatchSize * (Width - k * PatchSize); l++)
						{
							Data_of_Band[l] = kPatch.at<float>(0, l);

						}

						pBand->RasterIO(GF_Write, k * PatchSize, j * PatchSize, (Width - k * PatchSize), PatchSize, Data_of_Band,
							(Width - k * PatchSize), PatchSize, GDT_Float32, 0, 0);//闂備浇娉曢崰鎰板几婵犳艾绠?柣鎴ｅГ閺呮悂鏌￠崒妯衡偓鏍?偓姘鳖劑+1闂備浇娉曢崰鎰板几婵犳艾绠?柣鎴ｅГ閺呮悂鏌￠崒妯衡偓鏍?偓姘?秺閺屻劑鎮㈢拠娴嬫瀼缂佺偓婢橀ˇ鎵?偓姘?秺閺屻劑鎮㈤崨濠勪紕闂佸綊顥撻崗姗�寮?幘璇茬?闁汇値鍨辫ぐ褔鏌熸?搴″幋闁轰焦鎹囧?顒勫Χ閸℃?浼揇ata_of_Band

						delete[]Data_of_Band;
						_kPatch.release();
						kPatch.release();
					}
					else if (j == HLops - 1 && k != WLops - 1)
					{
						float* Data_of_Band = new float[(Height - j * PatchSize) * PatchSize];

						cv::Mat _kPatch;
						Image[i](cv::Rect(k* PatchSize, j* PatchSize, PatchSize, Height - j * PatchSize)).copyTo(_kPatch);
						cv::Mat kPatch = _kPatch.reshape(1, 1);

						for (int l = 0; l < (Height - j * PatchSize) * PatchSize; l++)
						{
							Data_of_Band[l] = kPatch.at<float>(0, l);

						}

						pBand->RasterIO(GF_Write, k * PatchSize, j * PatchSize, PatchSize, Height - j * PatchSize, Data_of_Band,
							PatchSize, Height - j * PatchSize, GDT_Float32, 0, 0);//闂備浇娉曢崰鎰板几婵犳艾绠?柣鎴ｅГ閺呮悂鏌￠崒妯衡偓鏍?偓姘鳖劑+1闂備浇娉曢崰鎰板几婵犳艾绠?柣鎴ｅГ閺呮悂鏌￠崒妯衡偓鏍?偓姘?秺閺屻劑鎮㈢拠娴嬫瀼缂佺偓婢橀ˇ鎵?偓姘?秺閺屻劑鎮㈤崨濠勪紕闂佸綊顥撻崗姗�寮?幘璇茬?闁汇値鍨辫ぐ褔鏌熸?搴″幋闁轰焦鎹囧?顒勫Χ閸℃?浼揇ata_of_Band

						delete[]Data_of_Band;
						_kPatch.release();
						kPatch.release();
					}
					else
					{
						float* Data_of_Band = new float[(Height - j * PatchSize) * (Width - k * PatchSize)];

						cv::Mat _kPatch;
						Image[i](cv::Rect(k* PatchSize, j* PatchSize, Width - k * PatchSize, Height - j * PatchSize)).copyTo(_kPatch);

						cv::Mat kPatch = _kPatch.reshape(1, 1);

						for (int l = 0; l < (Height - j * PatchSize) * (Width - k * PatchSize); l++)
						{
							Data_of_Band[l] = kPatch.at<float>(0, l);

						}

						pBand->RasterIO(GF_Write, k * PatchSize, j * PatchSize, Width - k * PatchSize, Height - j * PatchSize, Data_of_Band,
							Width - k * PatchSize, Height - j * PatchSize, GDT_Float32, 0, 0);//闂備浇娉曢崰鎰板几婵犳艾绠?柣鎴ｅГ閺呮悂鏌￠崒妯衡偓鏍?偓姘鳖劑+1闂備浇娉曢崰鎰板几婵犳艾绠?柣鎴ｅГ閺呮悂鏌￠崒妯衡偓鏍?偓姘?秺閺屻劑鎮㈢拠娴嬫瀼缂佺偓婢橀ˇ鎵?偓姘?秺閺屻劑鎮㈤崨濠勪紕闂佸綊顥撻崗姗�寮?幘璇茬?闁汇値鍨辫ぐ褔鏌熸?搴″幋闁轰焦鎹囧?顒勫Χ閸℃?浼揇ata_of_Band

						delete[]Data_of_Band;
						_kPatch.release();
						kPatch.release();
					}
				}
			}


		}
	}
	else
	{
		for (int i = 0; i < BandSize; i++)
		{
			GDALRasterBand* pBand = ImageData->GetRasterBand(i + 1);
			float* Data_of_Band = new float[Height * Width];

			cv::Mat strechimage = Image[i].reshape(1, 1);

			for (int j = 0; j < Height * Width; j++)
			{
				Data_of_Band[j] = strechimage.at<float>(0, j);

			}
			//闂備浇娉曢崰鏍?箺椤忓牆妫橀柕鍫濇?椤忓爼姊虹捄銊ユ瀺闁肩厧閰ｅΛ鍐?煛閸屾凹浼撻悗娈垮枛缁诲?鎷归悜顧奣_Float32闂備浇娉曢崰鎰板几婵犳艾绠?柧姘?�婚?閬嶆⒑鐠恒劌鏋戦柡瀣?煼楠炲繘鎮块?锝嗙煑闂備浇娉曢崰鎰板几婵犳艾绠?�瑰嫭澹嗙�癸拷
			pBand->RasterIO(GF_Write, 0, 0, Width, Height, Data_of_Band, Width, Height, GDT_Float32, 0, 0);
			delete[]Data_of_Band;

		}
	}




	
	delete[]Projection;

	//std::cout << "闂備浇娉曢崰鎰板几婵犳艾绠?柣鎴ｅГ閺呮悂鏌￠崒妯衡偓鏍?偓姘鳖攰閵囨劙鏌嗗?鍡樻?闂佸搫鍊堕崐鏍?偓姘?秺瀵?煡骞掑Δ浣规?闂佽法鍣﹂幏锟?!" << std::endl;
	GDALClose(ImageData);
}


// void BasicFunction::SaveImageAsTiff(std::vector<cv::Mat> Image, std::string SaveFilePath, char* Proj, double* GeoTransform)
// {
// 	  std::cout<<"692\n";
// 	const int BandSize =Image.size();
// 	const int Height = Image[0].rows;
// 	const int Width = Image[0].cols;

// 	GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
// 	GDALDataset* ImageData;
// 	//新建的格式为GDT_Float32，一定要对应
// 	ImageData = poDriver->Create(SaveFilePath.data(), Width, Height, BandSize, GDT_Float32, NULL);

// 	const char* Projection = Proj;
// 	ImageData->SetGeoTransform(GeoTransform);
// 	ImageData->SetProjection(Projection);

// 	int PatchSize = 1000;


// 	if (Width > 10 * PatchSize && Height > 10 * PatchSize)
// 	{
// 		int WLops = Width / PatchSize + 1;
// 		int HLops = Height / PatchSize + 1;

// 		//cout << float(Width) / float(PatchSize) <<"  " << Width / PatchSize << endl;
// 		//刚好凑整的时候，不加1
// 		if (float(Width) / float(PatchSize) == Width / PatchSize)
// 		{
// 			WLops = WLops - 1;
// 		}
// 		if (float(Height) / float(PatchSize) ==  Height / PatchSize)
// 		{
// 			HLops = HLops - 1;
// 		}
// 		//cout << WLops << "  " << HLops << endl;
// 		for (int i = 0; i < BandSize; i++)
// 		{
// 			GDALRasterBand* pBand = ImageData->GetRasterBand(i + 1);

// 			for (int j = 0; j < HLops; j++)
// 			{
// 				for (int k = 0; k < WLops; k++)
// 				{
// 					if (j != HLops - 1 && k != WLops - 1)
// 					{
// 						float* Data_of_Band = new float[PatchSize * PatchSize];

// 						cv::Mat _kPatch;
// 						Image[i](cv::Rect(k* PatchSize, j* PatchSize, PatchSize, PatchSize)).copyTo(_kPatch );
						
// 						cv::Mat kPatch = _kPatch.reshape(1, 1);

// 						for (int l = 0; l < PatchSize * PatchSize; l++)
// 						{
// 							Data_of_Band[l] = kPatch.at<float>(0, l);

// 						}
// 						pBand->RasterIO(GF_Write, k * PatchSize, j * PatchSize, PatchSize, PatchSize, Data_of_Band,
// 							PatchSize, PatchSize, GDT_Float32, 0, 0);//将第i+1个波段的数据存入Data_of_Band

// 						delete[]Data_of_Band;
// 						_kPatch.release();
// 						kPatch.release();
// 					}
// 					else if (j != HLops - 1 && k == WLops - 1)
// 					{
// 						float* Data_of_Band = new float[PatchSize * (Width - k * PatchSize)];

// 						cv::Mat _kPatch;
// 						Image[i](cv::Rect(k* PatchSize, j* PatchSize, (Width - k * PatchSize), PatchSize)).copyTo(_kPatch);

// 						cv::Mat kPatch = _kPatch.reshape(1, 1);
// 						for (int l = 0; l < PatchSize * (Width - k * PatchSize); l++)
// 						{
// 							Data_of_Band[l] = kPatch.at<float>(0, l);

// 						}

// 						pBand->RasterIO(GF_Write, k * PatchSize, j * PatchSize, (Width - k * PatchSize), PatchSize, Data_of_Band,
// 							(Width - k * PatchSize), PatchSize, GDT_Float32, 0, 0);//将第i+1个波段的数据存入Data_of_Band

// 						delete[]Data_of_Band;
// 						_kPatch.release();
// 						kPatch.release();
// 					}
// 					else if (j == HLops - 1 && k != WLops - 1)
// 					{
// 						float* Data_of_Band = new float[(Height - j * PatchSize) * PatchSize];

// 						cv::Mat _kPatch;
// 						Image[i](cv::Rect(k* PatchSize, j* PatchSize, PatchSize, Height - j * PatchSize)).copyTo(_kPatch);
// 						cv::Mat kPatch = _kPatch.reshape(1, 1);

// 						for (int l = 0; l < (Height - j * PatchSize) * PatchSize; l++)
// 						{
// 							Data_of_Band[l] = kPatch.at<float>(0, l);

// 						}

// 						pBand->RasterIO(GF_Write, k * PatchSize, j * PatchSize, PatchSize, Height - j * PatchSize, Data_of_Band,
// 							PatchSize, Height - j * PatchSize, GDT_Float32, 0, 0);//将第i+1个波段的数据存入Data_of_Band

// 						delete[]Data_of_Band;
// 						_kPatch.release();
// 						kPatch.release();
// 					}
// 					else
// 					{
// 						float* Data_of_Band = new float[(Height - j * PatchSize) * (Width - k * PatchSize)];

// 						cv::Mat _kPatch;
// 						Image[i](cv::Rect(k* PatchSize, j* PatchSize, Width - k * PatchSize, Height - j * PatchSize)).copyTo(_kPatch);

// 						cv::Mat kPatch = _kPatch.reshape(1, 1);

// 						for (int l = 0; l < (Height - j * PatchSize) * (Width - k * PatchSize); l++)
// 						{
// 							Data_of_Band[l] = kPatch.at<float>(0, l);

// 						}

// 						pBand->RasterIO(GF_Write, k * PatchSize, j * PatchSize, Width - k * PatchSize, Height - j * PatchSize, Data_of_Band,
// 							Width - k * PatchSize, Height - j * PatchSize, GDT_Float32, 0, 0);//将第i+1个波段的数据存入Data_of_Band

// 						delete[]Data_of_Band;
// 						_kPatch.release();
// 						kPatch.release();
// 					}
// 				}
// 			}


// 		}
// 	}
// 	else
// 	{
// 		for (int i = 0; i < BandSize; i++)
// 		{
// 			GDALRasterBand* pBand = ImageData->GetRasterBand(i + 1);
// 			float* Data_of_Band = new float[Height * Width];

// 			cv::Mat strechimage = Image[i].reshape(1, 1);

// 			for (int j = 0; j < Height * Width; j++)
// 			{
// 				Data_of_Band[j] = strechimage.at<float>(0, j);

// 			}
// 			//新建的格式为GDT_Float32，一定要对应
// 			pBand->RasterIO(GF_Write, 0, 0, Width, Height, Data_of_Band, Width, Height, GDT_Float32, 0, 0);
// 			delete[]Data_of_Band;

// 		}
// 	}




	
// 	delete[]Projection;

// 	//std::cout << "保存影像成功!" << std::endl;
// 	GDALClose(ImageData);
// }


double BasicFunction::CalEuDistance2D(double x1, double y1, double x2, double y2)
{
	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}


//根据多项式变化参数计算点的位置
cv::Point2f BasicFunction::PolynomialTransform(float x, float y, cv::Mat PolynomialParams)
{
	cv::Point2f Output;
	Output.x = PolynomialParams.at<float>(0, 0) + PolynomialParams.at<float>(1, 0) * x + PolynomialParams.at<float>(2, 0) * y +
		PolynomialParams.at<float>(3, 0) * x * x + PolynomialParams.at<float>(4, 0) * x * y + PolynomialParams.at<float>(5, 0) * y * y;
	Output.y = PolynomialParams.at<float>(6, 0) + PolynomialParams.at<float>(7, 0) * x + PolynomialParams.at<float>(8, 0) * y +
		PolynomialParams.at<float>(9, 0) * x * x + PolynomialParams.at<float>(10, 0) * x * y + PolynomialParams.at<float>(11, 0) * y * y;
	return Output;
}
//返回多项式校正的边界
DEMBoundary BasicFunction::GetBoundaryFromPolynomialParams(cv::Size ImageSize, cv::Mat PolynomialParams)
{
	Boundary ImageBoundary;

	cv::Point2f tmp;
	tmp = PolynomialTransform(0, 0, PolynomialParams);
	ImageBoundary.TopLeftLongitude = tmp.x;
	ImageBoundary.TopLeftLatitude = tmp.y;

	tmp = PolynomialTransform(float(ImageSize.width), 0, PolynomialParams);
	ImageBoundary.TopRightLongitude = tmp.x;
	ImageBoundary.TopRightLatitude = tmp.y;

	tmp = PolynomialTransform(0, float(ImageSize.height), PolynomialParams);
	ImageBoundary.BottomLeftLongitude = tmp.x;
	ImageBoundary.BottomLeftLatitude = tmp.y;

	tmp = PolynomialTransform(float(ImageSize.width), float(ImageSize.height), PolynomialParams);
	ImageBoundary.BottomRightLongitude = tmp.x;
	ImageBoundary.BottomRightLatitude = tmp.y;

	return BasicFunction::GetDEMBoundary(ImageBoundary);
}
//写入为指定格式的文件，除开jpeg2000
void BasicFunction::SaveImage(std::vector<cv::Mat> Image, std::string Format, std::string SaveFilePath, char* Proj, double* GeoTransform)
{
	const int BandSize = Image.size();
	const int Height = Image[0].rows;
	const int Width = Image[0].cols;

	//GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName(Format.c_str());

	GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName(Format.c_str());
	GDALDataset* ImageData;
	//新建的格式为GDT_Float32，一定要对应
	ImageData = poDriver->Create(SaveFilePath.data(), Width, Height, BandSize, GDT_Float32, NULL);

	const char* Projection = Proj;

	if (GeoTransform != nullptr)
	{
		ImageData->SetGeoTransform(GeoTransform);
	}
	if (Projection != nullptr)
	{
		ImageData->SetProjection(Projection);
	}

	for (int i = 0; i < BandSize; i++)
	{
		GDALRasterBand* poBand = ImageData->GetRasterBand(i + 1);
		float* Banddata = new float[Height * Width];

		cv::Mat strechimage = Image[i].reshape(1, 1);

		for (int j = 0; j < Height * Width; j++)
		{
			Banddata[j] = strechimage.at<float>(0, j);

		}
		//新建的格式为GDT_Float32，一定要对应
		poBand->RasterIO(GF_Write, 0, 0, Width, Height, Banddata, Width, Height, GDT_Float32, 0, 0);
		delete[]Banddata;

	}

	delete[]Projection;

	//std::cout << "保存影像成功!" << std::endl;
	GDALClose(ImageData);
}

void BasicFunction::SaveImageASOtherFormat(std::vector<cv::Mat> Image, std::string Format, std::string SaveFilePath, char* Proj, double* GeoTransform)
{
	GDALAllRegister();
	const int BandSize = Image.size();
	const int Height = Image[0].rows;
	const int Width = Image[0].cols;

	GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("MEM");
	GDALDataset* ImageData;
	//新建的格式为GDT_Float32，一定要对应
	ImageData = poDriver->Create("", Width, Height, BandSize, GDT_Int16, NULL);

	const char* Projection = Proj;
	if (GeoTransform != nullptr)
	{
		ImageData->SetGeoTransform(GeoTransform);
	}
	if (Projection != nullptr)
	{
		ImageData->SetProjection(Projection);
	}
	for (int i = 0; i < BandSize; i++)
	{
		GDALRasterBand* poBand = ImageData->GetRasterBand(i + 1);
		int* Banddata = new int[Height * Width];

		cv::Mat strechimage = Image[i].reshape(1, 1);

		for (int j = 0; j < Height * Width; j++)
		{
			Banddata[j] = strechimage.at<float>(0, j);

		}
		//新建的格式为GDT_Float32，一定要对应
		poBand->RasterIO(GF_Write, 0, 0, Width, Height, Banddata, Width, Height, GDT_Int16, 0, 0);
		delete[]Banddata;

	}
	GDALDriver* poDriverJP2 = GetGDALDriverManager()->GetDriverByName(Format.c_str());
	if (poDriverJP2 != nullptr)
	{
		poDriverJP2->CreateCopy(SaveFilePath.data(), ImageData, TRUE, 0, 0, 0);
	}
	else
	{
		std::cout << "没这个格式啊!" << std::endl;
	}
	//delete[]Projection;
	//GDALClose(ImageData);
	//delete[]Projection;

	//GDALClose(ImageData);
}
//GDAL格式转CV::Mat
bool BasicFunction::GDAL2Mat(GDALDataset* InputImageData, std::vector<cv::Mat>& OutputImage, char* Proj, double* GeoTransform)
{
	//获取投影转换
	const char* const_prj = InputImageData->GetProjectionRef();
	string tmp(const_prj);
	//转化一下
	ConstChar2Char(const_prj, Proj);


	//std::cout << Proj << std::endl;
	int Width = InputImageData->GetRasterXSize(); // 列 
	int Height = InputImageData->GetRasterYSize(); // 行 
	int BandSize = InputImageData->GetRasterCount();//波段数 

	double* adfGeoTransform = new double[6];
	InputImageData->GetGeoTransform(adfGeoTransform);

	GeoTransform = adfGeoTransform;
	for (int i = 0; i < 6; i++)
	{
		std::cout << GeoTransform[i] << std::endl;
	}

	std::vector <cv::Mat>Image;  // 定义元素为Mat的vecoter向量，向量的每一个元素存储一个波段的数据

	int PatchSize = 1000;


	if (Width > 10 * PatchSize && Height > 10 * PatchSize)
	{
		int WLops = Width / PatchSize;
		int HLops = Height / PatchSize;
		//刚好凑整的时候，不加1
		if (float(Width) / float(PatchSize) == Width / PatchSize)
		{
			WLops = WLops - 1;
		}
		if (float(Height) / float(PatchSize) == Height / PatchSize)
		{
			HLops = HLops - 1;
		}
		for (int i = 0; i < BandSize; i++)
		{
			GDALRasterBand* pBand = InputImageData->GetRasterBand(i + 1);//读取第i+1个波段的数据
			
			GDALDataType DataType = pBand->GetRasterDataType();

			cv::Mat Band_Image = cv::Mat(Height, Width, CV_32FC1);//将第i+1个波段的数据存存入A中

			for (int j = 0; j < HLops; j++)
			{
				for (int k = 0; k < WLops; k++)
				{
					if (j != HLops - 1 && k != WLops - 1)
					{
						float* Data_of_Band = new float[PatchSize * PatchSize];

						pBand->RasterIO(GF_Read, k * PatchSize, j * PatchSize, PatchSize, PatchSize, Data_of_Band,
							PatchSize, PatchSize, GDT_Float32, 0, 0);//将第i+1个波段的数据存入Data_of_Band

						
						cv::Mat kPatch(PatchSize, PatchSize, CV_32FC1, Data_of_Band);
						kPatch.copyTo(Band_Image(cv::Rect(k * PatchSize, j * PatchSize, PatchSize, PatchSize)));

						delete[]Data_of_Band;
					}
					else if (j != HLops - 1 && k == WLops - 1)
					{
						float* Data_of_Band = new float[PatchSize * (Width - k * PatchSize)];

						pBand->RasterIO(GF_Read, k * PatchSize, j * PatchSize, (Width - k * PatchSize), PatchSize, Data_of_Band,
							(Width - k * PatchSize), PatchSize, GDT_Float32, 0, 0);//将第i+1个波段的数据存入Data_of_Band
						cv::Mat kPatch(PatchSize, (Width - k * PatchSize), CV_32FC1, Data_of_Band);
						kPatch.copyTo(Band_Image(cv::Rect(k * PatchSize, j * PatchSize, (Width - k * PatchSize), PatchSize)));

						delete[]Data_of_Band;
					}
					else if (j == HLops - 1 && k != WLops - 1)
					{
						float* Data_of_Band = new float[(Height - j * PatchSize) * PatchSize];

						pBand->RasterIO(GF_Read, k * PatchSize, j * PatchSize, PatchSize, Height - j * PatchSize, Data_of_Band,
							PatchSize, Height - j * PatchSize, GDT_Float32, 0, 0);//将第i+1个波段的数据存入Data_of_Band
						cv::Mat kPatch(Height - j * PatchSize, PatchSize, CV_32FC1, Data_of_Band);
						kPatch.copyTo(Band_Image(cv::Rect(k * PatchSize, j * PatchSize, PatchSize, Height - j * PatchSize)));

						delete[]Data_of_Band;
					}
					else
					{
						float* Data_of_Band = new float[(Height - j * PatchSize) * (Width - k * PatchSize)];

						pBand->RasterIO(GF_Read, k * PatchSize, j * PatchSize, Width - k * PatchSize, Height - j * PatchSize, Data_of_Band,
							Width - k * PatchSize, Height - j * PatchSize, GDT_Float32, 0, 0);//将第i+1个波段的数据存入Data_of_Band
						cv::Mat kPatch(Height - j * PatchSize, Width - k * PatchSize, CV_32FC1, Data_of_Band);
						kPatch.copyTo(Band_Image(cv::Rect(k * PatchSize, j * PatchSize, Width - k * PatchSize, Height - j * PatchSize)));

						delete[]Data_of_Band;
					}
				}
			}

			Image.push_back(Band_Image.clone());
			if (i == BandSize - 1)
			{
				delete pBand;
			}
			Band_Image.release();
		}
	}
	else
	{
		for (int i = 0; i < BandSize; i++)
		{
			float* Data_of_Band = new float[Width * Height];

			GDALRasterBand* pBand = InputImageData->GetRasterBand(i + 1);//读取第i+1个波段的数据
			GDALDataType DataType = pBand->GetRasterDataType();

			pBand->RasterIO(GF_Read, 0, 0, Width, Height, Data_of_Band,
				Width, Height, GDT_Float32, 0, 0);//将第i+1个波段的数据存入pafScan
			cv::Mat Band_Image = cv::Mat(Height, Width, CV_32FC1, Data_of_Band);//将第i+1个波段的数据存存入A中
			Image.push_back(Band_Image.clone());

			if (i == BandSize - 1)
			{
				delete pBand;
			}
			delete[]Data_of_Band;
			Band_Image.release();
		}
	}
	OutputImage = Image;
	// GDALClose((GDALDatasetH)ImageData);
	return true;
}
//const char *转char *
void BasicFunction::ConstChar2Char(const char* Input, char* Output)
{
	/*string tmp(Input);
	Output = new char(tmp.length());
	strcpy(Output, tmp.c_str());
	std::cout << tmp << endl;
	cout << Output << endl;*/
	char* str = nullptr;
	str = const_cast<char*>(Input);
	int length = std::strlen(Input);
	for (int i = 0; i < length; i++) 
	{
		char itmp = str[i];
		Output[i] = itmp;
	}
}
//读取xml文件
bool BasicFunction::LoadXmlFile(std::string XmlFilePath,GF2Metadata & MetaData)
{


	tinyxml2::XMLDocument XmlFile;
	if (XmlFile.LoadFile(XmlFilePath.c_str()) != tinyxml2::XML_SUCCESS)
	{
		std::cout << "打开xml文件失败，请检查文件目录!" << std::endl;
		return false;
	}
	
	//读取根节点
	tinyxml2::XMLElement* rootElement = XmlFile.RootElement();
	const char* tMetadata;
	//波段数
	if (rootElement->FirstChildElement("Bands")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("Bands")->GetText();
		MetaData.Bands = atoi(tMetadata);
	}
	else
	{
		std::cout << "读取波段数失败!" << std::endl;
	}
	//左下角纬度
	if (rootElement->FirstChildElement("BottomLeftLatitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("BottomLeftLatitude")->GetText();
		MetaData.BottomLeftLatitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取左下角纬度失败!" << std::endl;
	}
	
	//左下角经度
	if (rootElement->FirstChildElement("BottomLeftLongitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("BottomLeftLongitude")->GetText();
		MetaData.BottomLeftLongitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取左下角经度失败!" << std::endl;
	}
	
	//右下角纬度
	if (rootElement->FirstChildElement("BottomRightLatitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("BottomRightLatitude")->GetText();
		MetaData.BottomRightLatitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取右下角纬度失败!" << std::endl;
	}
	
	//右下角经度
	if (rootElement->FirstChildElement("BottomRightLongitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("BottomRightLongitude")->GetText();
		MetaData.BottomRightLongitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取右下角经度失败!" << std::endl;
	}

	
	//椭球框架
	if (rootElement->FirstChildElement("EarthEllipsoid")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("EarthEllipsoid")->GetText();
		MetaData.EarthEllipsoid = tMetadata;
	}
	else
	{
		std::cout << "读取椭球框架失败!" << std::endl;
	}
	
	//地面分辨率
	if (rootElement->FirstChildElement("ImageGSD")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("ImageGSD")->GetText();
		MetaData.ImageGSD = atof(tMetadata);
	}
	else
	{
		std::cout << "读取地面分辨率失败!" << std::endl;
	}
	
	//等级
	if (rootElement->FirstChildElement("IntegrationLevel")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("IntegrationLevel")->GetText();
		MetaData.IntegrationLevel = tMetadata;
	}
	else
	{
		std::cout << "读取IntegrationLevel失败!" << std::endl;
	}
	
	//时间
	if (rootElement->FirstChildElement("IntegrationTime")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("IntegrationTime")->GetText();
		MetaData.IntegrationTime = atof(tMetadata);
	}
	else
	{
		std::cout << "读取IntegrationTime失败!" << std::endl;
	}
	
	//
	if (rootElement->FirstChildElement("MtfCorrection")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("MtfCorrection")->GetText();
		MetaData.MtfCorrection = tMetadata;
	}
	else
	{
		std::cout << "读取MtfCorrection失败!" << std::endl;
	}
	
	//
	if (rootElement->FirstChildElement("PitchSatelliteAngle")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("PitchSatelliteAngle")->GetText();
		MetaData.PitchSatelliteAngle = atof(tMetadata);
	}
	else
	{
		std::cout << "读取PitchSatelliteAngle失败!" << std::endl;
	}
	
	//
	if (rootElement->FirstChildElement("PitchViewingAngle")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("PitchViewingAngle")->GetText();
		MetaData.PitchViewingAngle = atof(tMetadata);
	}
	else
	{
		std::cout << "读取PitchViewingAngle失败!" << std::endl;
	}
	
	//格式
	if (rootElement->FirstChildElement("ProductFormat")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("ProductFormat")->GetText();
		MetaData.ProductFormat = tMetadata;
	}
	else
	{
		std::cout << "读取ProductFormat失败!" << std::endl;
	}

	
	//接收时间
	if (rootElement->FirstChildElement("ReceiveTime")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("ReceiveTime")->GetText();
		MetaData.ReceiveTime = tMetadata;
	}
	else
	{
		std::cout << "读取ReceiveTime失败!" << std::endl;
	}
	
	//
	if (rootElement->FirstChildElement("RollSatelliteAngle")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("RollSatelliteAngle")->GetText();
		MetaData.RollSatelliteAngle = atof(tMetadata);
	}
	else
	{
		std::cout << "读取RollSatelliteAngle失败!" << std::endl;
	}
	
	//
	if (rootElement->FirstChildElement("RollViewingAngle")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("RollViewingAngle")->GetText();
		MetaData.RollViewingAngle = atof(tMetadata);
	}
	else
	{
		std::cout << "读取RollViewingAngle失败!" << std::endl;
	}

	
	//
	if (rootElement->FirstChildElement("SatelliteAzimuth")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("SatelliteAzimuth")->GetText();
		MetaData.SatelliteAzimuth = atof(tMetadata);
	}
	else
	{
		std::cout << "读取SatelliteAzimuth失败!" << std::endl;
	}
	
	//卫星ID
	if (rootElement->FirstChildElement("SatelliteID")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("SatelliteID")->GetText();
		MetaData.SatelliteID = tMetadata;
	}
	else
	{
		std::cout << "读取卫星ID失败!" << std::endl;
	}

	//
	if (rootElement->FirstChildElement("SatelliteZenith")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("SatelliteZenith")->GetText();
		MetaData.SatelliteZenith = atof(tMetadata);
	}
	else
	{
		std::cout << "读取SatelliteZenith失败!" << std::endl;
	}
	
	//
	if (rootElement->FirstChildElement("SensorID")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("SensorID")->GetText();
		MetaData.SensorID = tMetadata;
	}
	else
	{
		std::cout << "读取SensorID失败!" << std::endl;
	}

	//
	if (rootElement->FirstChildElement("SolarAzimuth")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("SolarAzimuth")->GetText();
		MetaData.SolarAzimuth = atof(tMetadata);
	}
	else
	{
		std::cout << "读取SolarAzimuth失败!" << std::endl;
	}
	
	//
	if (rootElement->FirstChildElement("SolarZenith")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("SolarZenith")->GetText();
		MetaData.SolarZenith = atof(tMetadata);
	}
	else
	{
		std::cout << "读取SolarZenith失败!" << std::endl;
	}
	
	//左上角纬度
	if (rootElement->FirstChildElement("TopLeftLatitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("TopLeftLatitude")->GetText();
		MetaData.TopLeftLatitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取左上角纬度失败!" << std::endl;
	}
	
	//左上角经度
	if (rootElement->FirstChildElement("TopLeftLongitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("TopLeftLongitude")->GetText();
		MetaData.TopLeftLongitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取左上角经度失败!" << std::endl;
	}
	
	//右上角纬度
	if (rootElement->FirstChildElement("TopRightLatitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("TopRightLatitude")->GetText();
		MetaData.TopRightLatitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取右上角纬度失败!" << std::endl;
	}
	
	//右上角经度
	if (rootElement->FirstChildElement("TopRightLongitude")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("TopRightLongitude")->GetText();
		MetaData.TopRightLongitude = atof(tMetadata);
	}
	else
	{
		std::cout << "读取右上角经度失败!" << std::endl;
	}

	
	//
	if (rootElement->FirstChildElement("YawSatelliteAngle")->GetText() != NULL)
	{
		tMetadata = rootElement->FirstChildElement("YawSatelliteAngle")->GetText();
		MetaData.YawSatelliteAngle = atof(tMetadata);
	}
	else
	{
		std::cout << "读取YawSatelliteAngle失败!" << std::endl;
	}
	
	return 1;
}
//从元数据中获取影像边界
Boundary BasicFunction::GetBoundaryFromMetadata(GF2Metadata Metadata)
{
	Boundary OutputBoundary;
	OutputBoundary.BottomLeftLatitude = Metadata.BottomLeftLatitude;
	OutputBoundary.BottomLeftLongitude = Metadata.BottomLeftLongitude;
	OutputBoundary.BottomRightLatitude = Metadata.BottomRightLatitude;
	OutputBoundary.BottomRightLongitude = Metadata.BottomRightLongitude;
	OutputBoundary.TopLeftLatitude = Metadata.TopLeftLatitude;
	OutputBoundary.TopLeftLongitude = Metadata.TopLeftLongitude;
	OutputBoundary.TopRightLatitude = Metadata.TopRightLatitude;
	OutputBoundary.TopRightLongitude = Metadata.TopRightLongitude;
	return OutputBoundary;
}
//从GeotransformParams获取边界
Boundary BasicFunction::GetBoundaryFromGeotransformParams(double* GeotransformParams,cv::Size ImageSize)
{
	Boundary OutputBoundary;
	OutputBoundary.BottomLeftLatitude = GeotransformParams[3] + (0.5) * GeotransformParams[4] + (ImageSize.height - 0.5) * GeotransformParams[5];
	OutputBoundary.BottomLeftLongitude = GeotransformParams[0] + (0.5) * GeotransformParams[1] + (ImageSize.height - 0.5) * GeotransformParams[2];

	OutputBoundary.BottomRightLatitude = GeotransformParams[3] + (ImageSize.width - 0.5) * GeotransformParams[4] + (ImageSize.height - 0.5) * GeotransformParams[5];
	OutputBoundary.BottomRightLongitude = GeotransformParams[0] + (ImageSize.width - 0.5) * GeotransformParams[1] + (ImageSize.height - 0.5) * GeotransformParams[2];

	OutputBoundary.TopLeftLatitude = GeotransformParams[3] + (0.5) * GeotransformParams[4] + (0.5) * GeotransformParams[5];
	OutputBoundary.TopLeftLongitude = GeotransformParams[0] + (0.5) * GeotransformParams[1] + (0.5) * GeotransformParams[2];

	OutputBoundary.TopRightLatitude = GeotransformParams[3] + (ImageSize.width - 0.5) * GeotransformParams[4] + (0.5) * GeotransformParams[5];
	OutputBoundary.TopRightLongitude = GeotransformParams[0] + (ImageSize.width - 0.5) * GeotransformParams[1] + (0.5) * GeotransformParams[2];
	return OutputBoundary;
}
std::string BasicFunction::GetXmlFilePath(std::string ImageFilePath)
{
	//设置xml文件的路径
	std::string sxmlPath;
	sxmlPath.assign(ImageFilePath.begin(), (--ImageFilePath.end() - 4));
	sxmlPath = sxmlPath + ".xml";
	return sxmlPath;
}

//直方图匹配--------------------------------------------------------------------------------------
//寻找vector最小值位置
int BasicFunction::GetMinNumIndex(vector<float> data)
{
	float minnum(99999999);
	int index(0);
	for (int i(0); i < data.size(); ++i)
	{
		if (data[i] < minnum)
		{
			minnum = data[i];
			index = i;
		}
	}
	return index;
}

//cv::Mat BasicFunction::CHistMatch(cv::Mat img1, cv::Mat img2, int code)
//{
//	/// @brief img2到img1，单通道
//	/// @param img1 第一张影像
//	/// @param img2 第二张影像
//	/// @param code 数据类型----8(bit),12(bit),16(bit)
//	/// @return 
//	cv::Mat mimg1 = img1.clone();
//	cv::Mat mimg2 = img2.clone();
//	//变换数据类型
//	mimg1.convertTo(mimg1, CV_32FC1);
//	mimg2.convertTo(mimg2, CV_32FC1);
//	//变形
//	mimg1 = mimg1.reshape(1, mimg1.rows * mimg1.cols);
//	mimg2 = mimg2.reshape(1, mimg2.rows * mimg2.cols);
//	//灰度级
//	int length(1);
//	for (size_t i = 0; i < code; i++)
//	{
//		length *= 2;
//	}
//	//灰度直方图
//	std::vector<float> hist1(length), hist2(length);
//	long long int fucklaijisit = mimg1.rows;
//	for (long long int i = 0; i < fucklaijisit; i++)
//	{
//		int num1 = (int)mimg1.at<float>(i, 0);
//		if (num1 < 0)
//			num1 = 0;
//		if (num1 >= length)
//			num1 = length - 1;
//		int num2 = (int)mimg2.at<float>(i, 0);
//		if (num2 < 0)
//			num2 = 0;
//		if (num2 >= length)
//			num2 = length - 1;
//		hist1[num1] ++;
//		hist2[num2] ++;
//	}
//	for (size_t i = 0; i < length; i++)
//	{
//		hist1[i] /= mimg1.rows;
//		hist2[i] /= mimg2.rows;
//		//std::cout << hist1[i] << "  " << hist2[i] << std::endl;
//	}
//	//累计灰度直方图
//	std::vector<float> count_hist1 = hist1;
//	std::vector<float> count_hist2 = hist2;
//	for (size_t i = 1; i < count_hist1.size(); i++)
//	{
//		count_hist1[i] = count_hist1[i - 1] + count_hist1[i];
//	}
//	for (size_t i = 1; i < count_hist2.size(); i++)
//	{
//		count_hist2[i] = count_hist2[i - 1] + count_hist2[i];
//	}
//	//计算映射
//	vector<int> nindex;
//	for (int i(0); i < count_hist2.size(); ++i)
//	{
//		float num = count_hist2[i];
//		vector<float> absnum(count_hist1.size());
//		for (int j(0); j < count_hist1.size(); ++j)
//		{
//			absnum[j] = abs(count_hist1[j] - count_hist2[i]);
//		}
//		int index = GetMinNumIndex(absnum);
//		nindex.push_back(index);
//	}
//	//vector<int> nindex;
//	//for (int i(0); i < count_hist2.size(); ++i)
//	//{
//	//	float num = count_hist2[i];
//	//	int index(0);
//	//	for (int j(0); j < count_hist1.size()-1; ++j)
//	//	{
//	//		//两个值
//	//		float num1 = count_hist1[j];
//	//		float num2 = count_hist2[j + 1];
//	//		if (num <= num1)
//	//		{
//	//			index = j;
//	//		}
//	//		else if (num1 < num && num < num2)
//	//		{
//	//			float dert1 = abs(num-num1);
//	//			float dert2 = abs(num-num2);
//	//			index = dert1 > dert2 ? j + 1 : j;
//	//		}
//	//		else if (num >= num2)
//	//		{
//	//			index = j + 1;
//	//		}
//	//	}
//	//	nindex.push_back(index);
//	//}
//	//开始映射
//	for (long long int i = 0; i < fucklaijisit; i++)
//	{
//		int num = (int)mimg2.at<float>(i, 0);
//		if (num < 0)
//			num = 0;
//		if (num >= length)
//			num = length - 1;
//		mimg2.at<float>(i, 0) = nindex[num];
//	}
//	cv::Mat nimg2 = img2.clone();
//	mimg2.reshape(1, img2.rows).convertTo(nimg2, img2.type());
//	return nimg2;
//}

cv::Mat BasicFunction::CHistMatch(const cv::Mat img1, const cv::Mat img2)
{
	/// @brief img2到img1，单通道
	/// @param img1 第一张影像
	/// @param img2 第二张影像
	/// @return 
	cv::Mat mimg1 = img1.clone();
	cv::Mat mimg2 = img2.clone();
	//变换数据类型
	mimg1.convertTo(mimg1, CV_32FC1);
	mimg2.convertTo(mimg2, CV_32FC1);
	//变形
	mimg1 = mimg1.reshape(1, mimg1.rows * mimg1.cols);
	mimg2 = mimg2.reshape(1, mimg2.rows * mimg2.cols);
	//排序并记下序号
	vector<FloatSort> mmimg1(mimg1.rows), mmimg2(mimg2.rows);
	for (size_t i = 0; i < mmimg1.size(); i++)
	{
		mmimg1[i].index = i;
		mmimg1[i].value = mimg1.at<float>(i, 0);
		mmimg2[i].index = i;
		mmimg2[i].value = mimg2.at<float>(i, 0);
	}
	sort(mmimg1.begin(), mmimg1.end(), compare);
	sort(mmimg2.begin(), mmimg2.end(), compare);
	//根据排序结果进行映射
	for (size_t i = 0; i < mmimg1.size(); i++)
	{
		float num1 = mmimg1[i].value;
		float num2 = mmimg2[i].value;
		if (num1 != num2)
		{
			mimg2.at<float>(mmimg2[i].index, 0) = num1;
		}
	}
	//变换并返回
	cv::Mat nimg2 = cv::Mat(img2.rows, img2.cols, img2.type());
	mimg2.reshape(1, img2.rows).convertTo(nimg2, img2.type());
	return nimg2;
}

bool BasicFunction::compare(FloatSort a, FloatSort b)
{
	return a.value < b.value;
}

//把映射结果算一个直方图玩一玩
//int length = mmimg1[mmimg1.size() - 1].value >= mmimg2[mmimg2.size() - 1].value ? mmimg1[mmimg1.size() - 1].value : mmimg2[mmimg2.size() - 1].value;
//std::vector<float> hist1(length), hist2(length);
//long long int fucklaijisit = mimg1.rows;
//for (long long int i = 0; i < fucklaijisit; i++)
//{
//	int num1 = (int)mimg1.at<float>(i, 0);
//	if (num1 < 0)
//		num1 = 0;
//	if (num1 >= length)
//		num1 = length - 1;
//	int num2 = (int)mimg2.at<float>(i, 0);
//	if (num2 < 0)
//		num2 = 0;
//	if (num2 >= length)
//		num2 = length - 1;
//	hist1[num1] ++;
//	hist2[num2] ++;
//}
////保存
//ofstream fout1("testline1.txt");
//for (size_t i = 0; i < hist1.size(); i++)
//{
//	fout1 << hist1[i] << std::endl;
//}
//fout1.close();
//ofstream fout2("testline2.txt");
//for (size_t i = 0; i < hist2.size(); i++)
//{
//	fout2 << hist2[i] << std::endl;
//}
//fout2.close();

//template <typename T>
//vector<T> sort_indexes(vector<T>& v)
//{
//	vector<T> idx(v.size());
//	iota(idx.begin(), idx.end(), 0);
//	sort(idx.begin(), idx.end(),
//		[&v](T i1, T i2) {return v[i1] < v[i2]; });
//	return idx;
//}


cv::Mat BasicFunction::LinearStrech(cv::Mat img, float ratio, int code, int dtype)
{
	cv::Mat mimg = img.clone();
	mimg.convertTo(mimg, CV_32FC1);
	cv::Mat nimg = mimg.clone();
	int counts = img.rows * img.cols;
	mimg = mimg.reshape(1, counts);
	cv::sort(mimg, mimg, cv::SORT_EVERY_COLUMN + cv::SORT_ASCENDING);
	float cutmin = mimg.at<float>(int(counts * ratio), 0);
	float cutmax = mimg.at<float>(int(counts * (1 - ratio)) - 1, 0);
	int num(1);
	for (size_t i = 0; i < code; i++)
	{
		num *= 2;
	}
	nimg.setTo(cutmin, nimg < cutmin);
	nimg.setTo(cutmax, nimg > cutmax);
	nimg = (num - 1) * (nimg - cutmin) / (cutmax - cutmin);
	nimg.convertTo(nimg, dtype);
	return nimg;
}

cv::Point3f BasicFunction::WGS842Proj(int des_code, double lon, double lat,double H)
{
	std::cout<<"call WGS842Proj"<<std::endl;
	int src_code = 4326;
	std::cout<<"lon\t"<<lon<<"\tlat\t"<<lat<<std::endl;
	//转换-----------------------
	OGRSpatialReference source_geo, des_geo;
	OGRCoordinateTransformation* source2des;
	std::cout<<"des_code\t"<<des_code<<std::endl;
	source_geo.importFromEPSG(src_code);
	des_geo.importFromEPSG(des_code);
	char* src_coor = NULL;
	source_geo.exportToWkt(&src_coor);
	char* des_coor = NULL;
	des_geo.exportToWkt(&des_coor);
	source2des = OGRCreateCoordinateTransformation(&source_geo, &des_geo);
	std::cout<<"src_code\t"<<src_code<<std::endl;
	double x = lon;
	double y = lat;
	std::cout<<"x\t"<<x<<"\ty\t"<<y<<std::endl;
	source2des->Transform(1, &x, &y);

	std::cout<<"okkkkk"<<std::endl;
	return cv::Point3f(x, y, H);
}



int BasicFunction::search_epsg_code(string datum, string proj_mode, double lon)
{
	int code;
	// std::vector<std::string> baga;
	// Stringsplit(mode, "_", baga);
	// std::string datum = baga[0];
	// std::string proj_mode = baga[1];
	if (datum == "WGS84" and proj_mode == "6")
	{
		// WGS84 -> UTM projection
		if (lon >= 132 && lon <= 138) { code = 32653; }
		if (lon >= 126 && lon <= 132) { code = 32652; }
		if (lon >= 120 && lon <= 126) { code = 32651; }
		if (lon >= 114 && lon <= 120) { code = 32650; }
		if (lon >= 108 && lon <= 114) { code = 32649; }
		if (lon >= 102 && lon <= 108) { code = 32648; }
		if (lon >= 96 && lon <= 102) { code = 32647; }
		if (lon >= 90 && lon <= 96) { code = 32646; }
		if (lon >= 84 && lon <= 90) { code = 32645; }
		if (lon >= 78 && lon <= 84) { code = 32644; }
		if (lon >= 72 && lon <= 78) { code = 32643; }
	}
	if (datum == "CGCS2000" and proj_mode == "6")
	{
		if (lon >= 132 && lon <= 138) { code = 32653; }
		if (lon >= 126 && lon <= 132) { code = 32652; }
		if (lon >= 120 && lon <= 126) { code = 32651; }
		if (lon >= 114 && lon <= 120) { code = 32650; }
		if (lon >= 108 && lon <= 114) { code = 32649; }
		if (lon >= 102 && lon <= 108) { code = 32648; }
		if (lon >= 96 && lon <= 102) { code = 32647; }
		if (lon >= 90 && lon <= 96) { code = 32646; }
		if (lon >= 84 && lon <= 90) { code = 32645; }
		if (lon >= 78 && lon <= 84) { code = 32644; }
		if (lon >= 72 && lon <= 78) { code = 32643; }
	}
	if (datum == "CGCS2000" and proj_mode == "3")
	{
		if (lon >= 133.5) { code = 4554; }
		if (lon >= 130.5 && lon <= 133.5) { code = 4553; }
		if (lon >= 127.5 && lon <= 130.5) { code = 4552; }
		if (lon >= 124.5 && lon <= 127.5) { code = 4551; }
		if (lon >= 121.5 && lon <= 124.5) { code = 4550; }
		if (lon >= 118.5 && lon <= 121.5) { code = 4549; }
		if (lon >= 115.5 && lon <= 118.5) { code = 4548; }
		if (lon >= 112.5 && lon <= 115.5) { code = 4547; }
		if (lon >= 109.5 && lon <= 112.5) { code = 4546; }
		if (lon >= 106.5 && lon <= 109.5) { code = 4545; }
		if (lon >= 103.5 && lon <= 106.5) { code = 4544; }
		if (lon >= 100.5 && lon <= 103.5) { code = 4543; }
		if (lon >= 97.5 && lon <= 100.5) { code = 4542; }
		if (lon >= 94.5 && lon <= 97.5) { code = 4541; }
		if (lon >= 91.5 && lon <= 94.5) { code = 4540; }
		if (lon >= 88.5 && lon <= 91.5) { code = 4539; }
		if (lon >= 85.5 && lon <= 88.5) { code = 4538; }
		if (lon >= 82.5 && lon <= 85.5) { code = 4537; }
		if (lon >= 79.5 && lon <= 82.5) { code = 4536; }
		if (lon >= 76.5 && lon <= 79.5) { code = 4535; }
		if (lon <= 76.5) { code = 4534; }
	}
	if (datum == "CGCS2000" and proj_mode == "6")
	{
		if (lon >= 132) { code = 4512; }
		if (lon >= 126 && lon <= 132) { code = 4511; }
		if (lon >= 120 && lon <= 126) { code = 4510; }
		if (lon >= 114 && lon <= 120) { code = 4509; }
		if (lon >= 108 && lon <= 114) { code = 4508; }
		if (lon >= 102 && lon <= 108) { code = 4507; }
		if (lon >= 96 && lon <= 102) { code = 4506; }
		if (lon >= 90 && lon <= 96) { code = 4505; }
		if (lon >= 84 && lon <= 90) { code = 4504; }
		if (lon >= 78 && lon <= 84) { code = 4503; }
		if (lon <= 78) { code = 4502; }
	}
	if (datum == "Xian1980" and proj_mode == "3")
	{
		if (lon >= 133.5) { code = 2390; }
		if (lon >= 130.5 && lon <= 133.5) { code = 2389; }
		if (lon >= 127.5 && lon <= 130.5) { code = 2388; }
		if (lon >= 124.5 && lon <= 127.5) { code = 2387; }
		if (lon >= 121.5 && lon <= 124.5) { code = 2386; }
		if (lon >= 118.5 && lon <= 121.5) { code = 2385; }
		if (lon >= 115.5 && lon <= 118.5) { code = 2384; }
		if (lon >= 112.5 && lon <= 115.5) { code = 2383; }
		if (lon >= 109.5 && lon <= 112.5) { code = 2382; }
		if (lon >= 106.5 && lon <= 109.5) { code = 2381; }
		if (lon >= 103.5 && lon <= 106.5) { code = 2380; }
		if (lon >= 100.5 && lon <= 103.5) { code = 2379; }
		if (lon >= 97.5 && lon <= 100.5) { code = 2378; }
		if (lon >= 94.5 && lon <= 97.5) { code = 2377; }
		if (lon >= 91.5 && lon <= 94.5) { code = 2376; }
		if (lon >= 88.5 && lon <= 91.5) { code = 2375; }
		if (lon >= 85.5 && lon <= 88.5) { code = 2374; }
		if (lon >= 82.5 && lon <= 85.5) { code = 2373; }
		if (lon >= 79.5 && lon <= 82.5) { code = 2372; }
		if (lon >= 76.5 && lon <= 79.5) { code = 2371; }
		if (lon <= 76.5) { code = 2370; }
	}
	if (datum == "Xian1980" and proj_mode == "6")
	{
		if (lon >= 132) { code = 2348; }
		if (lon >= 126 && lon <= 132) { code = 2347; }
		if (lon >= 120 && lon <= 126) { code = 2346; }
		if (lon >= 114 && lon <= 120) { code = 2345; }
		if (lon >= 108 && lon <= 114) { code = 2344; }
		if (lon >= 102 && lon <= 108) { code = 2343; }
		if (lon >= 96 && lon <= 102) { code = 2342; }
		if (lon >= 90 && lon <= 96) { code = 2341; }
		if (lon >= 84 && lon <= 90) { code = 2340; }
		if (lon >= 78 && lon <= 84) { code = 2339; }
		if (lon <= 78) { code = 2338; }
	}
	if (datum == "Beijing1954" and proj_mode == "3")
	{
		if (lon >= 133.5) { code = 2442; }
		if (lon >= 130.5 && lon <= 133.5) { code = 2441; }
		if (lon >= 127.5 && lon <= 130.5) { code = 2440; }
		if (lon >= 124.5 && lon <= 127.5) { code = 2439; }
		if (lon >= 121.5 && lon <= 124.5) { code = 2438; }
		if (lon >= 118.5 && lon <= 121.5) { code = 2437; }
		if (lon >= 115.5 && lon <= 118.5) { code = 2436; }
		if (lon >= 112.5 && lon <= 115.5) { code = 2435; }
		if (lon >= 109.5 && lon <= 112.5) { code = 2434; }
		if (lon >= 106.5 && lon <= 109.5) { code = 2433; }
		if (lon >= 103.5 && lon <= 106.5) { code = 2432; }
		if (lon >= 100.5 && lon <= 103.5) { code = 2431; }
		if (lon >= 97.5 && lon <= 100.5) { code = 2430; }
		if (lon >= 94.5 && lon <= 97.5) { code = 2429; }
		if (lon >= 91.5 && lon <= 94.5) { code = 2428; }
		if (lon >= 88.5 && lon <= 91.5) { code = 2427; }
		if (lon >= 85.5 && lon <= 88.5) { code = 2426; }
		if (lon >= 82.5 && lon <= 85.5) { code = 2425; }
		if (lon >= 79.5 && lon <= 82.5) { code = 2424; }
		if (lon >= 76.5 && lon <= 79.5) { code = 2423; }
		if (lon <= 76.5) { code = 2422; }
	}
	if (datum == "Beijing1954" and proj_mode == "6")
	{
		if (lon >= 132) { code = 21463; }
		if (lon >= 126 && lon <= 132) { code = 21462; }
		if (lon >= 120 && lon <= 126) { code = 21461; }
		if (lon >= 114 && lon <= 120) { code = 21460; }
		if (lon >= 108 && lon <= 114) { code = 21459; }
		if (lon >= 102 && lon <= 108) { code = 21458; }
		if (lon >= 96 && lon <= 102) { code = 21457; }
		if (lon >= 90 && lon <= 96) { code = 21456; }
		if (lon >= 84 && lon <= 90) { code = 21455; }
		if (lon >= 78 && lon <= 84) { code = 21454; }
		if (lon <= 78) { code = 21453; }
	}
	return code;
}

void BasicFunction::Stringsplit(const std::string& str, const std::string& splits, std::vector<std::string>& res)
{
	if (str == "")		return;
	//在字符串末尾也加入分隔符，方便截取最后一段
	std::string strs = str + splits;
	size_t pos = strs.find(splits);
	int step = splits.size();
	res.clear();
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

std::vector<cv::Mat> BasicFunction::GetSubPatch(std::string Img_path, cv::Point2d TLPoint, cv::Point2d BRPoint)
{
	// 第一步：获取输入影像的仿射变换参数
	// 1、读取影像
	GDALRPCInfo RPCParams_ref;
	GDALDataset* InputImageData = BasicFunction::ReadGTiffImage(Img_path, RPCParams_ref);
	// 2、获取仿射变换参数
	double* adfGeoTransform = new double[6];
	InputImageData->GetGeoTransform(adfGeoTransform);

	// 第二步：根据给定的左上、右下角点，计算其在输入影像上的相方坐标
	// 1、获取坐标
	double TLP_x, TLP_y, BRP_x, BRP_y;
	TLP_x = TLPoint.x;
	TLP_y = TLPoint.y;
	BRP_x = BRPoint.x;
	BRP_y = BRPoint.y;

	// 2、根据仿射变换参数计算对应相方坐标，加减1是为了避免在int的操作中损失结果
	int Pix_x_l, Pix_x_r, Pix_y_l, Pix_y_r;
	Pix_x_l = int((TLP_x - adfGeoTransform[0]) / adfGeoTransform[1]) - 1;
	Pix_x_r = int((BRP_x - adfGeoTransform[0]) / adfGeoTransform[1]) + 1;
	Pix_y_l = int((TLP_y - adfGeoTransform[3]) / adfGeoTransform[5]) - 1;
	Pix_y_r = int((BRP_y - adfGeoTransform[3]) / adfGeoTransform[5]) + 1;

	int xsize = Pix_x_r - Pix_x_l;
	int ysize = Pix_y_r - Pix_y_l;

	// 第三步：根据相方坐标从输入影像中裁剪子影像
	int Width = InputImageData->GetRasterXSize(); // 列 
	int Height = InputImageData->GetRasterYSize(); // 行 
	int BandSize = InputImageData->GetRasterCount();//波段数

	std::vector<cv::Mat> sub_patch;
	for (int i = 0; i < BandSize; i++)
	{
		GDALRasterBand* pBand = InputImageData->GetRasterBand(i + 1);
		unsigned int* Data_of_Band = new unsigned int[xsize * ysize];
		pBand->RasterIO(GF_Read, Pix_x_l, Pix_y_l, xsize, ysize, Data_of_Band,
			xsize, ysize, GDT_Byte, 0, 0);//将第i+1个波段的数据存入Data_of_Band
		cv::Mat sub_band(ysize, xsize, CV_8UC1, Data_of_Band);

		sub_patch.push_back(sub_band.clone());
		delete[] Data_of_Band;
		if (i == BandSize - 1)
		{
			delete pBand;
		}
		sub_band.release();
	}
	//地理信息
	cv::Mat xxxx = cv::Mat::zeros(ysize, xsize, CV_64FC1);
	cv::Mat yyyy = cv::Mat::zeros(ysize, xsize, CV_64FC1);
	for (size_t i = 0; i < ysize; i++)
	{
		for (size_t j = 0; j < xsize; j++)
		{
			xxxx.at<double>(i, j) = TLP_x + j * adfGeoTransform[1];
			yyyy.at<double>(i, j) = TLP_y + i * adfGeoTransform[5];
		}
	}
	sub_patch.push_back(xxxx);
	sub_patch.push_back(yyyy);
	return sub_patch;
}


std::vector<cv::Mat> BasicFunction::GetSubPatch(std::string Img_path, cv::Point2d TLPoint, int width, int height)
{
	// 第一步：获取输入影像的仿射变换参数
	// 1、读取影像
	GDALRPCInfo RPCParams_ref;
	GDALDataset* InputImageData = BasicFunction::ReadGTiffImage(Img_path, RPCParams_ref);
	// 2、获取仿射变换参数
	double* adfGeoTransform = new double[6];
	InputImageData->GetGeoTransform(adfGeoTransform);

	// 第二步：根据给定的左上、右下角点，计算其在输入影像上的相方坐标
	// 1、获取坐标
	double TLP_x, TLP_y;
	TLP_x = TLPoint.x;
	TLP_y = TLPoint.y;

	// 2、根据仿射变换参数计算对应相方坐标，加减1是为了避免在int的操作中损失结果
	int Pix_x_l, Pix_y_l;
	Pix_x_l = int((TLP_x - adfGeoTransform[0]) / adfGeoTransform[1]) - 1;
	Pix_y_l = int((TLP_y - adfGeoTransform[3]) / adfGeoTransform[5]) - 1;

	int xsize = width;
	int ysize = height;

	// 第三步：根据相方坐标从输入影像中裁剪子影像
	int Width = InputImageData->GetRasterXSize(); // 列 
	int Height = InputImageData->GetRasterYSize(); // 行 
	int BandSize = InputImageData->GetRasterCount();//波段数

	std::vector<cv::Mat> sub_patch;
	for (int i = 0; i < BandSize; i++)
	{
		GDALRasterBand* pBand = InputImageData->GetRasterBand(i + 1);
		unsigned int* Data_of_Band = new unsigned int[xsize * ysize];
		pBand->RasterIO(GF_Read, Pix_x_l, Pix_y_l, xsize, ysize, Data_of_Band,
			xsize, ysize, GDT_Byte, 0, 0);//将第i+1个波段的数据存入Data_of_Band
		cv::Mat sub_band(ysize, xsize, CV_8UC1, Data_of_Band);

		sub_patch.push_back(sub_band.clone());
		delete[] Data_of_Band;
		if (i == BandSize - 1)
		{
			delete pBand;
		}
		sub_band.release();
	}
	//地理信息
	cv::Mat xxxx = cv::Mat::zeros(ysize, xsize, CV_64FC1);
	cv::Mat yyyy = cv::Mat::zeros(ysize, xsize, CV_64FC1);
	for (size_t i = 0; i < ysize; i++)
	{
		for (size_t j = 0; j < xsize; j++)
		{
			xxxx.at<double>(i, j) = TLP_x + j * adfGeoTransform[1];
			yyyy.at<double>(i, j) = TLP_y + i * adfGeoTransform[5];
		}
	}
	sub_patch.push_back(xxxx);
	sub_patch.push_back(yyyy);
	return sub_patch;
}

std::string BasicFunction::GetBaseMapMode(std::string imgpath)
{
	char* Proj;
	double* GeoTransform;
	//初始化GDAL库注册表
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	//读取影像
	GDALDataset* ImageData = (GDALDataset*)GDALOpen(imgpath.data(), GA_ReadOnly);
	//获取投影转换
	const char* const_prj = ImageData->GetProjectionRef();
	Proj = const_cast<char*>(const_prj);
	std::string tstr = (std::string)Proj;
	std::vector<std::string> linestrs;
	Stringsplit(tstr, "\"", linestrs);
	return linestrs[1];
}


int BasicFunction::readXML4Boundary(string xmlPath, vector<vector<double>> &CornerPoints )
{
    xmlDocPtr doc;           //定义解析文档指针
    xmlNodePtr curNode;      //定义结点指针(你需要它为了在各个结点间移动)
    xmlChar *szKey;          //临时字符串变量
    // char *szDocName="/home/wj/Downloads/ZY302/zy302a_bwd_016497_003126_20190520112310_01_sec_0001_1905228482.xml";
    char *szDocName=(char *)xmlPath.data();
    
    std::cout.precision(12);




    doc = xmlReadFile(szDocName,"GB2312",XML_PARSE_RECOVER); //解析文件
    //检查解析文档是否成功，如果不成功，libxml将指一个注册的错误并停止。
    //一个常见错误是不适当的编码。XML标准文档除了用UTF-8或UTF-16外还可用其它编码保存。
    //如果文档是这样，libxml将自动地为你转换到UTF-8。更多关于XML编码信息包含在XML标准中.
    if (NULL == doc)
    {  
       cout<<"Xml not parsed successfully\n";    
       return -1;
    }
    curNode = xmlDocGetRootElement(doc); //确定文档根元素
    /*检查确认当前文档中包含内容*/
    if (NULL == curNode)
    {
       cout<<"empty xml\n";
       xmlFreeDoc(doc);
       return -1;
    }
    /*在这个例子中，我们需要确认文档是正确的类型。“root”是在这个示例中使用文档的根类型。*/
    if (xmlStrcmp(curNode->name, BAD_CAST "sensor_corrected_metadata"))
    {
       cout<<"Xml of the wrong type, root node != sensor_corrected_metadata";
       xmlFreeDoc(doc);
       return -1;
    }

   
 
    curNode = curNode->xmlChildrenNode;
    while (curNode != NULL) {
    if ((!xmlStrcmp(curNode->name, (const xmlChar *)"productInfo"))) {
        curNode = curNode->xmlChildrenNode;
        break;
    }
    curNode= curNode->next;
    }

    while (curNode != NULL) {
    if ((!xmlStrcmp(curNode->name, (const xmlChar *)"ProductGeographicRange"))) {
        curNode = curNode->xmlChildrenNode;
        break;
    }
    curNode= curNode->next;
    }

   
    CornerPoints.clear();
    vector<double> LonLat;
    
    while (curNode != NULL) 
    {
        if ((!xmlStrcmp(curNode->name, (const xmlChar *)"LeftTopPoint"))) 
        {
            cout<<curNode->name<<endl;
            xmlNodePtr PointNode=curNode->xmlChildrenNode;
            while (PointNode != NULL) {
            if ((!xmlStrcmp(PointNode->name, (const xmlChar *)"Longtitude"))) 
            {
                cout<<"\t"<<PointNode->name<<"\t";
                szKey = xmlNodeGetContent(PointNode);
                cout<<szKey<<endl;
                LonLat.push_back(stof((char*)szKey));

            }
            else if ((!xmlStrcmp(PointNode->name, (const xmlChar *)"Latitude")))
            {
                cout<<"\t"<<PointNode->name<<"\t";
                szKey = xmlNodeGetContent(PointNode);
                cout<<szKey<<endl;
                LonLat.push_back(stof((char*)szKey));
                CornerPoints.push_back(LonLat);
                LonLat.clear();
            }
            PointNode= PointNode->next;
            }
        }
        else if ((!xmlStrcmp(curNode->name, (const xmlChar *)"RightTopPoint"))) 
        {
            cout<<curNode->name<<endl;
            xmlNodePtr PointNode=curNode->xmlChildrenNode;
            while (PointNode != NULL) {
                if ((!xmlStrcmp(PointNode->name, (const xmlChar *)"Longtitude"))) 
                {
                    cout<<"\t"<<PointNode->name<<"\t";
                    szKey = xmlNodeGetContent(PointNode);
                    cout<<szKey<<endl;
                    LonLat.push_back(stof((char*)szKey));

                }
                else if ((!xmlStrcmp(PointNode->name, (const xmlChar *)"Latitude")))
                {
                    cout<<"\t"<<PointNode->name<<"\t";
                    szKey = xmlNodeGetContent(PointNode);
                    cout<<szKey<<endl;
                    LonLat.push_back(stof((char*)szKey));
                    CornerPoints.push_back(LonLat);
                    LonLat.clear();
                }
                
                PointNode= PointNode->next;
                }
            }
        else if ((!xmlStrcmp(curNode->name, (const xmlChar *)"RightBottomPoint"))) 
        {
            cout<<curNode->name<<endl;
            xmlNodePtr PointNode=curNode->xmlChildrenNode;
            while (PointNode != NULL) {
                if ((!xmlStrcmp(PointNode->name, (const xmlChar *)"Longtitude"))) 
                {
                    cout<<"\t"<<PointNode->name<<"\t";
                    szKey = xmlNodeGetContent(PointNode);
                    cout<<szKey<<endl;
                    LonLat.push_back(stof((char*)szKey));

                }
                else if ((!xmlStrcmp(PointNode->name, (const xmlChar *)"Latitude")))
                {
                    cout<<"\t"<<PointNode->name<<"\t";
                    szKey = xmlNodeGetContent(PointNode);
                    cout<<szKey<<endl;
                    LonLat.push_back(stof((char*)szKey));
                    CornerPoints.push_back(LonLat);
                    LonLat.clear();
                }
                PointNode= PointNode->next;
                }
            }
        else if ((!xmlStrcmp(curNode->name, (const xmlChar *)"LeftBottomPoint"))) 
        {
            cout<<curNode->name<<endl;
            xmlNodePtr PointNode=curNode->xmlChildrenNode;
            while (PointNode != NULL) {
                if ((!xmlStrcmp(PointNode->name, (const xmlChar *)"Longtitude"))) 
                {
                    cout<<"\t"<<PointNode->name<<"\t";
                    szKey = xmlNodeGetContent(PointNode);
                    cout<<szKey<<endl;
                    LonLat.push_back(stof((char*)szKey));

                }
                else if ((!xmlStrcmp(PointNode->name, (const xmlChar *)"Latitude")))
                {
                    cout<<"\t"<<PointNode->name<<"\t";
                    szKey = xmlNodeGetContent(PointNode);
                    cout<<szKey<<endl;
                    LonLat.push_back(stof((char*)szKey));
                    CornerPoints.push_back(LonLat);
                    LonLat.clear();
                }
                PointNode= PointNode->next;
                }
        }
    curNode= curNode->next;

    }

	   

    return 0;
}
