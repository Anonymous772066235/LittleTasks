#include <sys/wait.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "cmdline.h"

#include "Reconstruction.h"
#include "BasicFunction.h"
#include "CRPCBased3DReconstruction.h"


//以下三个函数与循环创建目录有关
int _sysmkdir(const std::string& dir)
{
    int ret = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (ret && errno == EEXIST)
    {
        printf("dir[%s] already exist.\n",dir.c_str());
    }
    else if (ret)
    {
        printf("create dir[%s] error: %d %s\n" ,dir.c_str(),ret ,strerror(errno));
        return -1;
    }
    else
    {
        printf("create dir[%s] success.\n", dir.c_str());
    }
    return 0;
}

std::string __getParentDir(const std::string& dir)
{
    std::string pdir = dir;
    if(pdir.length() < 1 || (pdir[0] != '/')){
        return "";
    }
    while(pdir.length() > 1 && (pdir[pdir.length() -1] == '/')) pdir = pdir.substr(0,pdir.length() -1);

    pdir = pdir.substr(0,pdir.find_last_of('/'));
    return pdir;
}

int _sysmkdirs(const std::string& dir)
{
    int ret = 0;
    if(dir.empty())
        return -1;
    std::string pdir;
    if((ret = _sysmkdir(dir)) == -1){
        pdir = __getParentDir(dir);
        if((ret = _sysmkdirs(pdir)) == 0){
            ret = _sysmkdirs(dir);
        }
    }
    return ret;
}





int main(int argc, char** argv)
{

    CmdLine cmd;
   
   
    std::string uSer;
    std::string tIme;   
    std::string imgListJson;
    std::string intMaxParallax;
    std::string intMinParallax; 
    std::string intRowGridSize;
    std::string intColGridSize;
    std::string floatProjErrorThreshold;
    std::string floatDenoisingPara;

    cmd.add(make_option('u', uSer, "input_usrname"));
	cmd.add(make_option('t', tIme, "input_time"));
	cmd.add(make_option('i', imgListJson, "input_imagelist_json_format"));
	cmd.add(make_option('P', intMaxParallax, "input_Max_Parallax"));
    cmd.add(make_option('p', intMinParallax, "input_Min_Parallax"));
	cmd.add(make_option('r', intRowGridSize, "input_Row_GridSize"));
	cmd.add(make_option('c', intColGridSize, "input_Col_GridSize"));
	cmd.add(make_option('e', floatProjErrorThreshold,"input_Proj_Error_Threshold"));
    cmd.add(make_option('d', floatDenoisingPara,"input_Denoising_Para;"));


	try
	{
		if (argc == 1) throw std::string("Invalid parameter.");
		cmd.process(argc, argv);
	}
	catch (const std::string& s)
	{
		std::cerr << "Usage: " << argv[0] << '\n'
            <<"[-u]--input_usrname\n"
            <<"[-t]--input_time\n"
            <<"[-i]--input_imagelist_json_format\n"
            <<"[-P]--input_Max_Parallax\n"
            <<"[-p]--input_Min_Parallax\n"
            <<"[-r]--input_Row_GridSize\n"
            <<"[-c]--input_Col_GridSize\n"
            <<"[-e]--input_Proj_Error_Threshold\n"
            <<"[-d]--input_Denoising_Para;\n"
			<< std::endl;
		std::cerr << s << std::endl;
		return EXIT_FAILURE;
	}


    std::vector<string>  imgPaths_tmp;
    BasicFunction::Stringsplit(imgListJson,"]",imgPaths_tmp);
    std::string imgPaths=imgPaths_tmp[0];
    BasicFunction::Stringsplit(imgPaths,"[",imgPaths_tmp);
    imgPaths=imgPaths_tmp[1];
    // 3.以","分割整个字符串为vector

    std::vector<string>  imgPathVec;
    BasicFunction::Stringsplit(imgPaths,",",imgPathVec);
    
    

    // 4.去除引号，输出看看有没有问题
    std::cout<<std::endl<<"imgPaths"<<std::endl;
    for(int i=0;i<imgPathVec.size();i++)
    {
       
        BasicFunction::Stringsplit(imgPathVec[i],".",imgPaths_tmp);
        std::string tmp="";
        for(int j=0;j<imgPaths_tmp.size()-1;j++)
        {
        	tmp+=imgPaths_tmp[j]+".";
        }
        imgPathVec[i]=tmp.substr(0,tmp.length()-1);   //到这里为止，影像路径还是本地的
        std::cout<<imgPathVec[i]<<std::endl;
    }


    
    if(access((imgPathVec[0]+".xml").c_str(), F_OK ) != -1 )
    {
        std::cout<<std::endl<<(imgPathVec[0]+".xml")<<" exist."<<std::endl;
    }
    else
    {
        std::cout<<"no xml file! please add xml file."<<std::endl;
        system("pause");

    }
    
  





    //写出的反馈信息
    //工程目录

    std::string projectPath="/home/data/"+uSer+"/Reconstruction3D/"+tIme;

    std::string _Info2 = "\"task path\":\"" + projectPath + "\"";

    //创建日志文件

    std::string RunningStateMsg=BasicFunction::Convert2Json_new(200,"3DReconstruction is running",0,"false");

    std::cout<<std::endl<<"create folder"<<std::endl;
    if(_sysmkdirs(projectPath)!=0)
    {
        std::cout<<"failed when create folder"<<std::endl;
        system("pause");
    };


    BasicFunction::Write2Text(projectPath+"/RunningState_Log.txt",RunningStateMsg);

    //运行功能

    int MaxParallax=std::stoi(intMaxParallax);
    int MinParallax=std::stoi(intMinParallax);
    int RowGridSize=std::stoi(intRowGridSize);
    int ColGridSize=std::stoi(intColGridSize);
    float ProjErrorThreshold=std::stof(floatProjErrorThreshold);
    float DenoisingPara=std::stof(floatDenoisingPara);


    bool ok=true;
    if(MaxParallax<0 || MaxParallax>128)
    {

        ok=false;
        std::cout<<"\nThe range of MaxParallax should be 0 to 128.\n";
    }
    if(MinParallax<-128 || MinParallax>0)
    {

        ok=false;
        std::cout<<"The range of MinParallax should be -128 to 0.\n";
    }
    if(RowGridSize<1 || RowGridSize>100 || ColGridSize<1 || ColGridSize>100 )
    {

        ok=false;
        std::cout<<"The range of RowGridSize and ColGridSize should be 1 to 100.\n";
    }
    if(ProjErrorThreshold<1 || ProjErrorThreshold>20 )
    {

        ok=false;
        std::cout<<"The range of ProjErrorThreshold should be 1 to 20.\n";
    }
    if(DenoisingPara<1 || DenoisingPara>10 )
    {

        ok=false;
        std::cout<<"The range of DenoisingPara should be 1 to 10.\n";
    }
    



    if(ok)
    {
        
        CRPCBased3DReconstruction func;
        func.SetParams(imgPathVec,MinParallax,MaxParallax,RowGridSize,ColGridSize,ProjErrorThreshold,DenoisingPara,projectPath);

        // std::string pointCloudSavePath=projectPath+"/pointCould.txt";
        func.run();

        RunningStateMsg=BasicFunction::Convert2Json_new(200,"3DReconstruction is done",1,"true");
        BasicFunction::Write2Text(projectPath+"/RunningState_Log.txt",RunningStateMsg);

        std::cout<<RunningStateMsg<<std::endl;
        const int MAX_LENGTH = 512;
        char buffer[MAX_LENGTH];
        getcwd(buffer, 512);
        std::string current_path = buffer;

        BasicFunction::Write2Text(current_path+"/PointCloud_Path.txt",projectPath+"/PointCloud.txt");

        
    }
    else{
        RunningStateMsg=BasicFunction::Convert2Json_new(500,"3DReconstruction failed",1,"true");
        BasicFunction::Write2Text(projectPath+"/RunningState_Log.txt",RunningStateMsg);
        std::cout<<RunningStateMsg<<std::endl;
    }	
	
    std::cout<<_Info2<<std::endl;
	
    
    return 0;

}