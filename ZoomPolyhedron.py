# File     :ZoomPolyhedron.py
# Author   :WooChi
# Time     :2021/08/24
# Function :1)对一个近似立方体(或其他凸多面体)缩放，返回新坐标
# Version  :1.0

import numpy as np
from scipy.spatial import ConvexHull


class Polyhedron():
    vertexs = np.array(
        [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]])
    volume = 0
    centroid = np.array([0, 0, 0])

    def __init__(self, vertexs):
        self.vertexs = vertexs
        self.volume = ConvexHull(vertexs).volume
        self.centroid = np.mean(self.vertexs, axis=0)

    def zoom(self, z=1):
        if z > 0:
            newVertexs = (self.vertexs - self.centroid) * pow(z, 1 / 3) + self.centroid
            return newVertexs
        else:
            print('缩放尺度应大于0!!!')

    def show(self):
        print('\n/-------------properties of polyhedron-------------->')
        print('-------- vertexs --------')
        print(self.vertexs)
        print('-------- volume --------')
        print(self.volume)
        print('-------- centroid --------')
        print(self.centroid)
        print('<--------------properties of polyhedron--------------/\n')


if __name__ == '__main__':
    # Cube01_vts = np.array(
    #     [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]])
    Cube01_vts = np.array(
        [[0, 0, 0], [0, 1.2, 0], [1, 6, 0], [1, 1, 0], [0, 0, 1], [2, 0, 1], [0, 1, 1.5], [1, 1, 1]])
    # 7个点也可
    # Cube01_vts = np.array(
    #     [[0, 0, 0], [0, 1.2, 0], [1, 6, 0], [1, 1, 0], [0, 0, 1], [2, 0, 1], [0, 1, 1.5]])
    # 6个点也可
    # Cube01_vts = np.array(
    #     [[0, 0, 0], [0, 1.2, 0], [1, 6, 0], [1, 1, 0], [0, 0, 1], [2, 0, 1]])

    zm = 1.1
    # zm = 0.5
    # zm = 35
    # zm = -0.5

    Cube01 = Polyhedron(Cube01_vts)
    Cube02_vts = Cube01.zoom(zm)
    Cube02 = Polyhedron(Cube02_vts)

    Cube01.show()
    Cube02.show()

    print(Cube02.volume / Cube01.volume)
