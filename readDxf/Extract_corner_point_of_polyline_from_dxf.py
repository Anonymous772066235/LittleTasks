# File     :Extract_corner_point_of_polyline_from_dxf.py
# Author   :WooChi
# Time     :2022/03/14
# Function :从dxf中提取多段线（房屋、可残缺）的角点
# Version  :

import time
import numpy as np
import matplotlib.pyplot as plt
from icecream import ic
import dxfgrabber


def getPointsFromDxf(path):
    """从def格式dlg中提取房屋拐点"""
    dxf = dxfgrabber.readfile(path)
    corner_point = []
    allpoints = []
    plt.figure(figsize=(9, 9),dpi=200)
    for e in dxf.entities:
        if e.dxftype == 'POLYLINE':
            print(e.dxftype, e.layer)
            points = e.points
            if e.is_closed:
                # 多段线闭合，
                allpoints.extend(points)

                angle = []
                # 计算第一个点所在的角
                vector1 = np.array(points[1]) - np.array(points[0])
                vector2 = np.array(points[0]) - np.array(points[-1])
                cos_angle = vector1.dot(vector2) / (
                        np.linalg.norm(vector1) * np.linalg.norm(vector2))
                angle.append(abs(cos_angle))
                # 计算中间点所在角
                for i in range(len(points) - 2):
                    vector1 = np.array(points[i + 1]) - np.array(points[i])
                    vector2 = np.array(points[i + 2]) - np.array(points[i + 1])
                    cos_angle = vector1.dot(vector2) / (
                            np.linalg.norm(vector1) * np.linalg.norm(vector2))
                    angle.append(abs(cos_angle))
                # 计算最后一点所在角
                vector1 = np.array(points[0]) - np.array(points[-1])
                vector2 = np.array(points[-1]) - np.array(points[-2])
                cos_angle = vector1.dot(vector2) / (
                        np.linalg.norm(vector1) * np.linalg.norm(vector2))
                angle.append(abs(cos_angle))
                # 根据角度判别点是否为拐点
                judge = np.array(angle) < 0.98
                for i in range(len(judge)):
                    if judge[i]:
                        corner_point.append(points[i])

                # 为了绘图，增加一点
                points.append(e.points[0])
            else:
                # 若多段线没有闭合，则只用计算中间点所在角，首尾去除
                allpoints.extend(points)

                angle = [1]
                for i in range(len(points) - 2):
                    vector1 = np.array(points[i + 1]) - np.array(points[i])
                    vector2 = np.array(points[i + 2]) - np.array(points[i + 1])
                    cos_angle = vector1.dot(vector2) / (
                            np.linalg.norm(vector1) * np.linalg.norm(vector2))
                    angle.append(abs(cos_angle))
                angle.append(1)

                judge = np.array(angle) < 0.98
                for i in range(len(judge)):
                    if judge[i]:
                        corner_point.append(points[i])

            points = np.array(points)
            # 绘制房屋轮廓
            plt.plot(points[:, 0], points[:, 1])

    allpoints = np.array(allpoints)
    ic(len(allpoints))
    ic(len(corner_point))
    corner_point = np.unique(corner_point, axis=0)
    ic(len(corner_point))
    # 绘制所有读取出的未经筛选的所有点
    plt.scatter(allpoints[:, 0], allpoints[:, 1], marker='o', color='b', label='all ' + str(len(allpoints)))
    # 绘制筛选的点
    plt.scatter(corner_point[:, 0], corner_point[:, 1], marker='^', color='r', label='all ' + str(len(corner_point)))
    plt.legend()
    plt.axis('equal')
    plt.pause(20)
    return corner_point[:, :2]


def run():
    t_start = time.time()
    # dxfpath
    # dxfpath = "D:/Program Files/JetBrains/PycharmProjects/TLS_DLG_COPY_20211227/data/PLUS/0302/S5_CUT.dxf"
    dxfpath = "D:/Program Files/JetBrains/PycharmProjects/TLS_DLG_COPY_20211227/data/PCB/PCB.dxf"
    # dxfpath = "D:/Program Files/JetBrains/PycharmProjects/TLS_DLG_COPY_20211227/data/NEW/C2.dxf"

    corner_point = getPointsFromDxf(dxfpath)
    ic(corner_point)
    # path2save = "D:/Program Files/JetBrains/PycharmProjects/Plane2Point/output/PLUS_DLG.txt"
    # path2save = "D:/Program Files/JetBrains/PycharmProjects/Plane2Point/output/PCB_DLG.txt"
    # path2save = "D:/Program Files/JetBrains/PycharmProjects/Plane2Point/output/NEW_DLG.txt"
    # np.savetxt(path2save,corner_point,delimiter=',')

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
