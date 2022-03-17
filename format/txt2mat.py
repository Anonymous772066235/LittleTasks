# File     :txt2mat.py
# Author   :WooChi
# Time     :2021/11/29
# Function :txt转mat
# Version  :

import os
import time
import numpy as np
import scipy.io as sio
from icecream import ic

def txt2mat(inpath):
    '''txt转mat（涉及列之间的交换、当前路径与工作路径转换、获取当前路径下文件列表并进行判断以及批处理）'''
    path_origin = os.getcwd()  # 当前路径
    os.chdir(inpath)  # 切换到工作路径
    for file in os.listdir():  # 读取工作路径下的文件列表
        if file[-3:] == "csv":  # 判断文件格式
            outname = file[:-4] + '.mat'  # 创建mat文件名
            # data = np.loadtxt(file, delimiter=',', usecols=(0, 1, 4))  # 按列进行数据读取
            data = np.loadtxt(file, delimiter=',')  # 按列进行数据读取
            # data[:, [0, 1]] = data[:, [1, 0]]  # 两列数据互换
            ic(data)
            sio.savemat(outname, {'data': data})  # 保存为mat文件
            print(file + '-->' + outname + '\tSuccessfully!')  # 日志输出
    os.chdir(path_origin)  # 切回当前路径


def run():
    '''输入+运行+计时'''
    t_start = time.time()

    inpath = 'D:/Program Files/JetBrains/PycharmProjects/BathymetricDepthModel/output/20200719'
    inpath = 'D:/Program Files/JetBrains/PycharmProjects/BathymetricDepthModel/data/ICESat-2/oahu/01/20190613'
    txt2mat(inpath)

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
