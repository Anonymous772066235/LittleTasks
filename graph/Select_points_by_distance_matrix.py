# File     :Select_points_by_distance_matrix.py
# Author   :WooChi
# Time     :2022/03/23
# Function :通过图的距离矩阵来寻求两个点集的重合点
# Version  :

import os
import math
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from icecream import ic
import Ransac_match as rsm


def DistanceMatrice(Pset):
    """计算距离矩阵"""
    L = np.zeros((len(Pset), len(Pset)))
    for i in range(0, len(Pset)):
        for j in range(0, len(Pset)):
            if i != j:
                L[i, j] = np.around(np.linalg.norm(Pset[i, :] - Pset[j, :]), decimals=3)
    for i in range(0, len(Pset)):
        L[i, i] = 0
    return L


def CompareDM(DM1, DM2, threshod=0.5):
    """比较距离矩阵各行之间的相似度"""
    # 相似度矩阵（m*n)
    M = np.zeros(shape=(len(DM2), len(DM1)))
    for i in range(len(DM2)):
        S = []
        for j in range(len(DM1)):
            score = 0
            begin = 0
            for m in range(len(DM2[0, :])):
                for n in range(begin, len(DM1[0, :])):
                    if abs(DM2[i][m] - DM1[j][n]) < threshod:
                        score += (threshod - abs(DM2[i][m] - DM1[j][n])) / threshod
                        begin = n + 1
                        break

            M[i, j] = score
            S.append(score)
    return M


def findRepeatNumber(nums):
    """求数组里的重复值"""
    rpt = []
    for i in range(len(nums) - 1):
        if nums[i] in nums[i + 1:]:
            rpt.append(nums[i])
    return np.unique(rpt)


def get_index(lst=None, item=''):
    """求数组里的某元素的所有索引"""
    return [i for i in range(len(lst)) if lst[i] == item]


def matchM(M):
    """根据相似度矩阵进行匹配"""
    argmax_M = np.argmax(M, axis=1)
    max_M = np.max(M, axis=1)
    ic(argmax_M)
    ic(max_M)
    repeat = findRepeatNumber(argmax_M)
    ic(repeat)
    for i in range(len(repeat)):
        todo = get_index(argmax_M, repeat[i])
        ic(todo)
        todomax = max_M[todo]
        ic(todomax)
        turemaxidx = todo[np.argmax(todomax)]
        ic(turemaxidx)
        for j in range(len(todo)):
            if todo[j] != turemaxidx:
                ic(todo[j])
                while not (np.argmax(M[todo[j], :]) in argmax_M):
                    M[todo[j], np.argmax(M[todo[j]])] = 0
                argmax_M[todo[j]] = np.argmax(M[todo[j]])
    return argmax_M


def VisualizeMacthLine(P_dlg, P_dopp, row_ind, col_ind):
    """画线"""
    for i in range(len(row_ind)):
        X = [P_dlg[row_ind, 0], P_dopp[col_ind, 0]]
        Y = [P_dlg[row_ind, 1], P_dopp[col_ind, 1]]
        plt.plot(X, Y, color='b')


def transformation(DOPP, R, T):
    import numpy as np
    dopp_t = np.matrix(np.dot(R, DOPP.transpose()) + T)
    dopp = dopp_t.transpose()
    return np.array(dopp)


def fun():
    points1 = np.loadtxt("D:/Program Files/JetBrains/PycharmProjects/Plane2Point/output/ok/NEW_DLG.txt", delimiter=',')
    points2 = np.loadtxt("D:/Program Files/JetBrains/PycharmProjects/Plane2Point/output/ok/intersection_new.txt",
                         delimiter=',')

    L1 = np.sort(DistanceMatrice(points1), axis=1, kind='quicksort', order=None)[:, 1:]
    L2 = np.sort(DistanceMatrice(points2), axis=1, kind='quicksort', order=None)[:, 1:]

    if len(points1) < len(points2):
        points1, points2 = points2, points1
        L1, L2 = L2, L1
    _, M = CompareDM(L1, L2, 1)
    ic(M)
    ic(len(M))

    argmax = matchM(M)
    ic(argmax)

    P = points1[argmax]

    R, T, rpt = rsm.ca_rt(points2, P, 8, 100, 3)
    ic(R)
    ic(T)
    points2 = rsm.transformation(points2, R, T)

    row = np.arange(len(points2))
    plt.figure(figsize=(9, 9), dpi=150)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    plt.scatter(points1[:, 0], points1[:, 1], marker='.', color='black', label='dlg ' + str(len(points1)))
    plt.scatter(points2[:, 0], points2[:, 1], color='r', label='pc ' + str(len(points2)))
    plt.scatter(P[:, 0], P[:, 1], color='c', label='dlg_select ' + str(len(P)))

    VisualizeMacthLine(points2, points1, row, argmax)
    plt.axis('equal')
    plt.xlabel('X', fontsize=22)
    plt.ylabel('Y', fontsize=22)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.title(str(rpt))
    plt.legend(loc='best')
    plt.show()

    np.savetxt("D:/Program Files/JetBrains/PycharmProjects/Plane2Point/output/ok/NEW_DLG_select.txt", P, delimiter=',')


def run():
    t_start = time.time()

    fun()

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
