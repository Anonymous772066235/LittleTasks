#include "get_cfg.h"
 
#include <fstream>
#include <iostream>
 
using namespace std;
/*desc:
这里通过map实现了cfg文件中kv对的存取，方便操作。
注意这里的文件读取的操作，已经getline等相关函数。
注意整体的逻辑： 打开文件、从上到下获取行、注释忽略、建立map、存取kv对，这样，我们就可以得到cfg文件中所有有用的参数了。
最后FindInConfig函数中，我们将找到的kv输出，当然，在使用过程中，我们可以直接获取进行进一步的计算工作。
*/
 
bool IsSpace(char c)
{
	if (c == ' ' || c == '\t')
	{
		return true;
	}
	else
	{
		return false;
	}
}
 
bool IsCommentChar(char c)
{
	if (c == COMMENT_CHAR)
	{
		return true;
	}
	else
	{
		return false;
	}
}
 
// trim函数的作用是把一个字符串左边和右边的空格去掉，即为trim
void Trim(string & str) // 引用传参，这样在函数中修改该参数也会修改相应的变量
{
	if (str.empty())
	{
		return;
	}
	int i, start_pos, end_pos;
	for (i = 0; i < str.size(); i++)
	{
		if (!IsSpace(str[i]))
		{
			break;//跳出
		}
	}
	if (i == str.size())//如果该行全是空格，则该行最后一个字符为"\n"，此时i == str.size()
	{
		str = "";
		return;
	}
	start_pos = i; // 获取到非空格的初始位置
 
	for (i = str.size() - 1; i >= 0; i--)
	{
		if (!IsSpace(str[i]))
		{
			break;
		}
	}
	end_pos = i;
	str = str.substr(start_pos, end_pos - start_pos + 1);
}
 
bool AnalyseLine(const string & line, string & key, string & value) // 分析一行，如果是注释行，则不处理，如果是k-v行，则提取出key-value值。
{
	if (line.empty())
	{
		return false;
	}
	int start_pos = 0, end_pos = line.size() - 1, pos;
	if ((pos = line.find(COMMENT_CHAR)) != -1)
	{
		if (0 == pos)
		{
			return false; // 如果一行的开头是#,说明是注释，则 不需要
		}
		end_pos = pos - 1; // 可能是注释在k-v后的情况
	}
	string new_line = line.substr(start_pos, end_pos - start_pos + 1); // 删掉后半部分的注释 FIX_ME： 这里应该是减错了吧
	// 下面pos的赋值时必要的，这样，就可在后面得到Key和value值了。
	if ((pos = new_line.find("=")) == -1) //说明前面没有 = 号
	{
		return false;
	}
	key = new_line.substr(0, pos); // 获得key
	value = new_line.substr(pos + 1, end_pos + 1 - (pos + 1)); // 获得value
	Trim(key);
	if (key.empty())
	{
		return false;
	}
	Trim(value); // 因为这里的key和value都是引用传递，可以直接被修改，所以不用返回
	return true;
}
 
 
// 读取一个cfg文件并保存到map中,map中的key和value都是string类型的。
bool ReadConfig(const string & filename, map<string, string> & m)
{
	m.clear(); // 删除map容器中的所有k-v对
	ifstream infile(filename.c_str());
	if (!infile)
	{
		cout << "file open failed!" << endl; // 文件打开失败，返回错误信息。
		return false;
	}
	string line, key, value; // 为了后续保存kv对
	while (getline(infile, line))
	{
		if (AnalyseLine(line, key, value))
		{
			m[key] = value; // 保存到map容器中的方法。
		}
	}
	infile.close(); // 当读取完之后，就要关闭这个文件。
	return true;
}
 
void PrintConfig(const map<string, string> & m)
{
	map<string, string>::const_iterator mite;
	for (mite = m.begin(); mite != m.end(); mite++)
	{
		cout << mite->first << "=" << mite->second << endl;
	}
}
 
void FindInConfig(map<string, string>  m, string  key) // 注意：之前用的一直都是string类型，所以这里用的也是string key,而不是char key。
{
	map<string, string>::iterator it;
	it = m.find(key);
	if (it == m.end())
	{
		cout << "there is no " << key << endl;
	}
	else
	{
		cout << it->second << endl;
	}
}
