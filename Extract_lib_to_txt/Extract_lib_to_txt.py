# File      :Extract_lib_to_txt.py
# Auther    :WooChi
# Time      :2022/06/08
# Version   :1.0
# Function  :自动提取lib文件夹下的.lib文件，根据文件名分别导出到release_lib.txt和debug_lib.txt中

import os

def Extract_lib(path):
    release = open("release_lib.txt", 'w')
    debug = open("debug_lib.txt", 'w')
    for file in os.listdir(path):
        if file[-5:]=='d.lib':
            debug.write(file)
            debug.write('\n')
        elif file[-4:]=='.lib':
            release.write(file)
            release.write('\n')
    debug.close()
    release.close()

if __name__ == '__main__':
    path=r"D:\Program Files (x86)\OPenCV\x64\vc16\lib"
    Extract_lib(path)
    pass
