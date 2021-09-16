# File      :demH2.py
# Author    :WJ
# Function  :
# Time      :2021/08/08
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
        if idx1[0] <0.001:
            h.append(target[i, 2] - data[idx1[1], 2])
    h = np.array(h)
    mean_h = np.mean(h)
    return mean_h, h


def heightDifference_1(target, data):
    h = target[:, 2] - data[:, 2]
    mean_h = np.mean(h)
    return mean_h, h


if __name__ == "__main__":
    start = time.time()
    # dem_20200719_20180704
    # dem_20190721_20180704
    # dem_20190524_20180619
    # dem_20190421_20190301
    # dem_20190222_20190224

    # ph=np.loadtxt('D:/Program Files/JetBrains/PycharmProjects/GUI_BathymetricModel/data/ATL03/20190421_order_validate/seafloor_validate.txt')
    # ph=ph[:,[0,1,4]]
    # ic(ph)

    dem01 = np.loadtxt(
        'D:/Program Files/JetBrains/PycharmProjects/GUI_BathymetricModel/Output/0823/dem_20190222_20190224_test.txt',
        delimiter=' ')
    dem02 = np.loadtxt(
        'D:/Program Files/JetBrains/PycharmProjects/GUI_BathymetricModel/Output/0823/dem_20190222_20190224_validate.txt',
        delimiter=' ')
    # ic(len(ph))
    ic(len(dem01))
    ic(len(dem02))

    if len(dem01)>=len(dem02):
        rdm = np.random.randint(0, len(dem02), 1000)
        # ic(rdm)
        mean, h = heightDifference_2(dem02[rdm], dem01,k=1)
    else:
        rdm = np.random.randint(0, len(dem01), 1000)
        # ic(rdm)
        mean, h = heightDifference_2(dem01[rdm], dem02, k=1)

    # mean,h = heightDifference_2(ph, dem01,k=1)

    ic(h)
    ic(mean)

    ic('Running time:', time.time() - start)
    # 验证（有没有加错，即数据上移还是下移）
