#include<iostream>
#include<vector>
#include<string>
//  sudo apt-get install libxml2 
//  sudo ln -s /usr/include/libxml2/libxml   /usr/include/libxml
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
using namespace std;

namespace xmlRead{



/// \param depth  递归深度
/// \param _xmlNodePtr 节点对象指针
void ParserXML(int &depth,xmlNodePtr _xmlNodePtr){
    bool ishas_child=false;
    // 所有子节点
    xmlNodePtr xmlNodePtr1= _xmlNodePtr->children;
    int count=0;
    while (xmlNodePtr1){
        // 判断节点类型
        if(xmlNodePtr1->type!=XML_TEXT_NODE){
        // xmlStrcmp(xmlNodePtr1->name,BAD_CAST "text")
            count++;
            // 子节点个数
            int childEleCount= xmlChildElementCount(xmlNodePtr1);
            cout.width(depth);
            if(childEleCount==0){
                ishas_child=false;
                // 表明无子节点
                cout<<count<<"."<<xmlNodePtr1->name<<" --> "<<xmlNodeGetContent(xmlNodePtr1)<<endl;
            }else{
                ishas_child=true;
                // 表明有子节点
                cout<<count<<"."<<xmlNodePtr1->name<<endl;
            }
            // 遍历节点属性
            xmlAttr* xmlAttr1=xmlNodePtr1->properties;
            if(xmlAttr1){
                cout.width(depth+1);
                cout<<"=>";
                while (xmlAttr1!=NULL){
                    // 判断使用存在属性 ， xmlGetProp获取属性值
                    if(xmlHasProp(xmlNodePtr1,xmlAttr1->name)){
                        cout<<" "<<xmlAttr1->name<<":"<<
                            xmlGetProp(xmlNodePtr1,xmlAttr1->name);
                    }
                    xmlAttr1=xmlAttr1->next;
                }
                cout<<endl;
            }
            // 递归调用
            if(ishas_child){
                depth+=10;
                ParserXML(depth,xmlNodePtr1);
            }
        }
        xmlNodePtr1=xmlNodePtr1->next;
    }
    if(depth>0){
        depth-=10;
    }
}

int  readXML(string xmlPath, vector<vector<double>> CornerPoints )
{
    xmlDocPtr doc;           //定义解析文档指针
    xmlNodePtr curNode;      //定义结点指针(你需要它为了在各个结点间移动)
    xmlChar *szKey;          //临时字符串变量
    // char *szDocName="/home/wj/Downloads/ZY302/zy302a_bwd_016497_003126_20190520112310_01_sec_0001_1905228482.xml";
    char *szDocName=(char *)xmlPath.data();
    
    cout.precision(12);




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

    int deep=0;
    ParserXML(deep,curNode);
    cout<<deep<<endl;

   
 
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

    for(int i=0;i<CornerPoints.size();i++)
    {
        
        cout<<CornerPoints[i][0]<<"\t"<<CornerPoints[i][1]<<endl;
        
    }
    

    return 0;
}

}
