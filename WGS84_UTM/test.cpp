#include"COORCONV.h"

#include<iostream>
using namespace std;


int main()
{
    WGS84Corr t;
    UTMCoor test;
    t.lat = 39.225697;
    t.log = 115.292040;
    int zone =floor(t.log/6)+31;
	LatLonToUTMXY(DegToRad(t.lat), DegToRad(t.log),zone,test);

    std::cout.precision(8);
	std::cout << test.x << "  " << test.y << std::endl;

    WGS84Corr t2;
	
	//长沙的zone为49  长沙的中央子午线经度为111，推出zone为49
	//return DegToRad(-183.0 + (zone * 6.0));
	UTMXYToLatLon(test.x, test.y, zone, false,t2);
	double lat = RadToDeg(t2.lat);  //注意弧度转角度！！！
	double log = RadToDeg(t2.log);
    std::cout << lat << "  " << log << std::endl;
    return 0;

}