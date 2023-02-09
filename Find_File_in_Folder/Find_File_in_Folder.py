# File      :Find_File_in_Folder.py
# Auther    :WooChi
# Time      :2023/02/01
# Version   :1.0
# Function  :

import os
import shutil

def Find_and_Move_File(srcPath,tarPath,fileSuffix):

    g = os.walk(srcPath)
    for path, dir_list, file_list in g:
        for file_name in file_list:
            if file_name[-len(fileSuffix):]==fileSuffix:
                print(os.path.join(path, file_name))
                shutil.move(os.path.join(path, file_name), os.path.join(tarPath, file_name))



if __name__ == '__main__':
    srcPath=r""
    tarPath=r""
    Find_and_Move_File(srcPath,tarPath,".mp4")
