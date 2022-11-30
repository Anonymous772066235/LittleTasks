#include "CRPCBased3DReconstruction.h"
#include <iostream>
//设置参数
void CRPCBased3DReconstruction::SetParams( std::vector<std::string>  ImageListName, float Mind, float Maxd, int Gridr, int Gridc, int ProjectErr_T, int m,std::string projpath)
{

	p_ImageListName = ImageListName;
	p_Mind = Mind;
	p_Maxd = Maxd;
	p_Gridr = Gridr;
	p_Gridc = Gridc;
	p_ProjectErr_T = ProjectErr_T;
	p_m = m;
	p_projpath=projpath;
}
//运行三维地形重建
void CRPCBased3DReconstruction::run()
{
	runSamePointTrackSGM();

	runMulti_SpaceIntersection();
}
//提取同名点
void CRPCBased3DReconstruction::
runSamePointTrackSGM()
{
	//初始化
	CSamePointTrackSGM func;

	//设置输入参数
	func.create(p_projpath,p_ImageListName, p_Mind, p_Maxd,  p_Gridr, p_Gridc);
	func.compute();
	CorrespondigPoints = func.GetSamePoints();
}
//多片前方交会
void CRPCBased3DReconstruction::runMulti_SpaceIntersection()
{
	int k = 20;
	vector<SATPoint3D> Points;
	CMulti_SpaceIntersection func;

	func.SetInputPathParams(p_ImageListName);
	
	func.SetInputCorreImagePoints(CorrespondigPoints);
	func.SetOtherParams(p_ProjectErr_T, k, p_m,p_projpath);

	std::string RunningStateMsg=BasicFunction::Convert2Json_new(200,"3DReconstruction is running",0.8,"false");
	BasicFunction::Write2Text(p_projpath+"/RunningState_Log.txt",RunningStateMsg);

	func.Multi_SpaceIntersection(Points,p_projpath+"/PointCloud.txt");

}
