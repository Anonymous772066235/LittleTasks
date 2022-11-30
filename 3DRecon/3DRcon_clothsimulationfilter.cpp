#include <sys/wait.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "cmdline.h"

#include "PointCloud2DEM.h"
#include "BasicFunction.h"



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
    std::string pointCloudPath;
    std::string floatCSF_Size;
    std::string floatClassificationThreshold;

    cmd.add(make_option('u', uSer, "input_usrname"));
	cmd.add(make_option('t', tIme, "input_time"));
	cmd.add(make_option('p', pointCloudPath, "input_pointCloud_Path"));
	cmd.add(make_option('g', floatCSF_Size, "input_CSF_Grid_Size"));
	cmd.add(make_option('c', floatClassificationThreshold, "input_Classification_Threshold"));


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
            <<"[-p]--input_pointCloud_Path\n"
            <<"[-g]--input_CSF_Grid_Size\n"
            <<"[-c]--input_Classification_Threshold\n"
			<< std::endl;
		std::cerr << s << std::endl;
		return EXIT_FAILURE;
	}


    //工程目录

    std::string projectPath="/home/data/"+uSer+"/Reconstruction3D/"+tIme;

    std::string _Info2 = "\"task path\":\"" + projectPath + "\"";

    //创建日志文件

    std::string RunningStateMsg=BasicFunction::Convert2Json_new(200,"CSF_Filter is running",0,"false");

    if(_sysmkdirs(projectPath)!=0)
    {
        std::cout<<"failed when create folder"<<std::endl;
        system("pause");
    };


    BasicFunction::Write2Text(projectPath+"/RunningState_Log.txt",RunningStateMsg);

    //运行功能




    float CSF_Size=std::stof(floatCSF_Size);
    float ClassificationThreshold=std::stof(floatClassificationThreshold);


    bool ok=true;
    if(CSF_Size<0.1 || CSF_Size>100 )
    {

        ok=false;
        std::cout<<"The range of CSF_Size should be 0.1 to 100.\n";
    }
    if(ClassificationThreshold<0.1 ||ClassificationThreshold>10 )
    {

        ok=false;
        std::cout<<"The range of ClassificationThreshold should be 0.1 to 10.\n";
    }


    if(ok)
    {
        PointCloud2DEM::Points2GroundPoints(pointCloudPath,projectPath,true, CSF_Size, ClassificationThreshold);
        RunningStateMsg=BasicFunction::Convert2Json_new(200,"CSF_Filter is done",1,"true");
        BasicFunction::Write2Text(projectPath+"/RunningState_Log.txt",RunningStateMsg);
        std::cout<<RunningStateMsg<<std::endl;
        const int MAX_LENGTH = 512;
        char buffer[MAX_LENGTH];
        getcwd(buffer, 512);
        std::string current_path = buffer;

        BasicFunction::Write2Text(current_path+"/GroundPointCloud_Path.txt",projectPath+"/GroundPointCloud.txt");
    }
    else{
        RunningStateMsg=BasicFunction::Convert2Json_new(500,"CSF_Filter failed",1,"true");
        BasicFunction::Write2Text(projectPath+"/RunningState_Log.txt",RunningStateMsg);
        std::cout<<RunningStateMsg<<std::endl;
    }

    

    std::cout<<_Info2<<std::endl;

    return 0;
}
