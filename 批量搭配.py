# File     :批量搭配.py
# Author   :WooChi
# Time     :2021/12/16
# Function :输入两目录A,B，分别获取两目录下的文件，AB相互搭配
# Version  :1.0

import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from icecream import ic


def BatchMatch(directory_0, directory_1):
    ic()
    file_list_0 = os.listdir(directory_0)
    for i in range(len(file_list_0)):
        file_list_0[i] = os.path.join(directory_0, file_list_0[i])

    file_list_1 = os.listdir(directory_1)
    for i in range(len(file_list_1)):
        file_list_1[i] = os.path.join(directory_1, file_list_1[i])



    Match = []
    for i in range(len(file_list_0)):
        for j in range(i, len(file_list_1)):
            Match.append([i, j])
    Match=np.array(Match)
    # print(file_list_0)
    # print(file_list_1)
    # print(Match)
    # for i in range(len(Match[:, 0])):
    #     print(Match[i, :])
    #     ic(file_list_0[Match[i, 0]], file_list_1[Match[i, 1]])
    return file_list_0, file_list_1, Match


def run():
    t_start = time.time()

    directory_0 = 'F:\Chrome Downloads\sentinel-2_oahu'
    directory_1 = 'F:\Chrome Downloads\Oahu'
    file_list_0, file_list_1, Match=BatchMatch(directory_0, directory_1)
    print(file_list_0)
    print(file_list_1)
    print(Match)
    for i in range(len(Match[:,0])):
        print(Match[i,:])
        ic(file_list_0[Match[i,0]],file_list_1[Match[i,1]])


    # <editor-fold desc="time output">
    t_consume = time.time() - t_start
    h = t_consume // 3600
    m = (t_consume - h * 3600) // 60
    s = t_consume - h * 3600 - m * 60
    print('---------------------------------------------------')
    print('Time consuming: %d hours %d minutes %.3f seconds' % (h, m, s))
    print('---------------------------------------------------')
    # </editor-fold>


if __name__ == '__main__':
    run()
