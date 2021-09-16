# File      :demH.py
# Author    :WJ
# Function  :
# Time      :2021/08/07
# Version   :
# Amend     :


import time
from scipy.spatial import KDTree
from icecream import ic
import sys
sys.setrecursionlimit(1000000)
import numpy as np
np.set_printoptions(precision=3, suppress=True, threshold=np.inf, linewidth=100)


def denoising(data, sigma=0.2):
    toDelete = []
    for i in range(len(data)):
        mean_temp = np.mean(data[i, :])
        if np.max(data[i, :]) - mean_temp > sigma or mean_temp - np.min(data[i, :]) > sigma:
            toDelete.append(i)
    data2 = np.delete(data, toDelete, axis=0)
    return data2


def heightDifference_2(target, data, k=10):
    """target:控制点点云  data:用以构建树结构的地面点云  k:平面方向上最近点的个数"""
    # 1.找最近点
    tree = KDTree(data[:, :2], 2)
    h = []
    r = 0.001
    for i in range(len(target)):
        idx1 = tree.query(target[i, :2], k=k, distance_upper_bound=r)
        ic(idx1)
        if idx1[0] == 0.001:
            h.append(target[i, 2] - data[idx1[1], 2])
    h = np.array(h)
    mean_h=np.mean(h)
    return mean_h,h


def heightDifference_1(target, data):
    h=target[:,2]-data[:,2]
    mean_h = np.mean(h)
    return mean_h, h


if __name__ == "__main__":
    start = time.time()


    dem01 = np.loadtxt('D:/Program Files (x86)/PyCharm/PycharmProjects/PythonProject/GUI_BathymetricModel/testdixing_20190222&20180316.txt', delimiter=' ')
    dem02 = np.loadtxt('D:/Program Files (x86)/PyCharm/PycharmProjects/PythonProject/GUI_BathymetricModel/testdixing_20190524&20180316.txt', delimiter=' ')


    
    ic(len(dem01))
    ic(len(dem02))
    rdm = np.random.randint(0,len(dem01), 100)
    mean,h = heightDifference_1(dem01[rdm], dem02[rdm])
    # mean,h = heightDifference_1(dem01, dem02)

    ic(h)
    ic(mean)

    ic('Running time:', time.time() - start)
    # 验证（有没有加错，即数据上移还是下移）
