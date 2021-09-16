# File      :area_poly.py
# Author    :WJ
# Function  :计算多变形集1与多边形集2之间的交集的面积
# Time      :2021/08/13
# Version   :
# Amend     :

from shapely.geometry import Polygon, Point
import matplotlib.pyplot as plt
import numpy as np


def intersection_2poly(poly_1, poly_2):
    '''计算两个多边形交集的面积'''
    poly_1 = Polygon(poly_1)
    poly_2 = Polygon(poly_2)
    in_s = poly_1.intersection(poly_2).area
    draw_2poly(poly_1, poly_2, in_s) # 可视化两个多边形
    return in_s


def intersection_polys(polys1, polys2):
    '''计算多边形集的交集的面积，返回二维数组，以polys1的多边形个数为行，polys2的个数为列'''
    area = np.zeros((len(polys1), len(polys2)))
    for i in range(len(polys1)):
        for j in range(len(polys2)):
            area[i, j] = intersection_2poly(polys1[i], polys2[j])
    return area


def draw_2poly(poly_1, poly_2, in_s):
    '''绘制出两个多边形'''
    poly_1 = Polygon(poly_1)
    poly_2 = Polygon(poly_2)

    edgeXY_T = np.array(poly_1.exterior.xy)
    edgeXY = edgeXY_T.T

    edgeXY_T2 = np.array(poly_2.exterior.xy)
    edgeXY2 = edgeXY_T2.T

    plt.plot(edgeXY[:, 0], edgeXY[:, 1])
    plt.plot(edgeXY2[:, 0], edgeXY2[:, 1])
    plt.title(str(in_s))
    plt.pause(5)
    plt.close()


if __name__ == '__main__':
    a = [[1, 2], [1, 3], [2, 4], [6, 3]]
    b = [[2, 1], [3, 3], [4, 4], [5, 1]]
    c = [[1, 1], [0, 0], [-1, 0]]

    d = [a]
    e = [a, b, c]

    # 输入均为多边形集(三维列表):[[多边形1],[多边形2]...[多边形n]]或[[多边形]]
    # 其中多边形(二维列表):[[点1],[点2],[点3],[点1]]
    # 其中点(二维坐标):[x,y]
    area = intersection_polys(d, e)
    print(area)
    print()
    # area = intersection_polys(e, e)
    # print(area)
