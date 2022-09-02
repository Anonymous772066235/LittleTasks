# File      :Construct_and_visualize_graph_structure_by_distance_matrix.py
# Auther    :WooChi
# Time      :2022/08/24
# Version   :1.0
# Function  :根据加权距离矩阵构图绘图，求节点间最短距离

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from icecream import ic
import matplotlib.patches as mpatches
from scipy.optimize import linear_sum_assignment

def readAnnex3(path):
    data3 = pd.read_csv(path, sep=',', header=None, encoding='unicode_escape', skiprows=2, index_col=0)
    return data3, data3.to_numpy()


def drawNet(Matrix, path):
    plt.figure(figsize=(12, 12))
    G = nx.Graph()  # 建立一个空的无向图G
    n = len(Matrix)

    for i in range(n):
        for j in range(i, n):
            G.add_nodes_from([i + 1], name=i + 1, kind=3, color='dodgerblue')
            if Matrix[i][j] > 0 and Matrix[i][j] != np.inf:
                G.add_edge(i + 1, j + 1, weight=Matrix[i][j])

    for i in range(n):
        # 1类用户
        if i in [17, 18, 19, 21, 27, 28]:
            attrs = {i: {"kind": 1, "color": 'coral'}}
            nx.set_node_attributes(G, attrs)
        # 2类用户
        if i in [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 22, 23, 24]:
            attrs = {i: {"kind": 2, "color": 'gold'}}
            nx.set_node_attributes(G, attrs)

        # 快速充电站
        if i in [6, 13, 18, 22]:
            attrs = {i: {"kind": 0, "color": 'lime'}}
            nx.set_node_attributes(G, attrs)

    # 移动储能车存放点
    attrs = {1: {"kind": np.inf, "color": 'red'}}
    nx.set_node_attributes(G, attrs)
    # 台区变
    attrs = {5: {"kind": -1, "color": 'gray'}, 14: {"kind": -2, "color": 'gray'}, 19: {"kind": -3, "color": 'gray'},
             27: {"kind": -4, "color": 'gray'}}
    nx.set_node_attributes(G, attrs)

    # 绘图
    nx.draw_networkx_edge_labels(G, pos=nx.kamada_kawai_layout(G), edge_labels=nx.get_edge_attributes(G, "weight"),
                                 font_size=8)
    nx.draw(G, pos=nx.kamada_kawai_layout(G), with_labels=True, font_weight='bold')
    nx.draw_networkx_nodes(G, pos=nx.kamada_kawai_layout(G),
                           node_color=nx.get_node_attributes(G, "color").values())

    plt.rcParams['font.sans-serif'] = ['SimHei']  # 中文字体设置-黑体
    plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像是负号'-'显示为方块的问题

    labels = ['存放点', '充电站', '台变区', '一类点', '二类点', '三类点']
    color = ["red", "lime", "gray", 'coral', 'gold', 'dodgerblue']



    patches = [mpatches.Patch(color=color[i], label="{:s}".format(labels[i])) for i in range(len(color))]
    ax = plt.gca()
    # 下面一行中bbox_to_anchor指定了legend的位置
    ax.legend(handles=patches, ncol=1, loc=3)  # 生成legend

    plt.savefig(path, dpi=300, bbox_inches='tight', pad_inches=0)
    # plt.pause(1)
    plt.close()
    return G


if __name__ == '__main__':
    path_annex3 = "附件3：29个节点数据.csv"
    data_pd, data_np = readAnnex3(path_annex3)

    path_netpic = "net.png"
    G=drawNet(data_np, path_netpic)



    TBQ=[5,14,19,27]
    CDZ=[6,13,18,22]
    # 路径
    X=[]
    # 距离
    S=[0,0,0,0]
    # 返程最短距离矩阵(需要进行二分匹配)
    s = np.zeros((4, 4))


    for n,i in enumerate(TBQ):
        # 先确定往程
        X.append(nx.dijkstra_path(G, 1, i))
        S[n]+=nx.dijkstra_path_length(G, 1, i)
        # 再确定往程，需要方案进行分配
        for m,j in  enumerate(CDZ):
            s[n,m]=nx.dijkstra_path_length(G,i, j)+nx.dijkstra_path_length(G,j, 1)

    ic(s)
    # 获得充电站分配方案
    row_ind, col_ind = linear_sum_assignment(s)

    ic(row_ind)
    ic(col_ind)
    # 根据分配方案确定每个台变区车辆去往的充电站
    for i in range(len(row_ind)):
        S[row_ind[i]]+=s[row_ind[i],col_ind[i]]
        X[row_ind[i]].extend(nx.dijkstra_path(G,TBQ[row_ind[i]], CDZ[col_ind[i]])[1:])
        X[row_ind[i]].extend(nx.dijkstra_path(G, CDZ[col_ind[i]], 1)[1:])

    ic(X)
    ic(S)

    # s: array([[7.9, 21.1, 27.3, 17.6],
    #           [13.6, 16., 27.5, 18.7],
    #           [19.2, 31.3, 16., 15.8],
    #           [17.8, 25.4, 18.1, 11.7]])
    # row_ind: array([0, 1, 2, 3])
    # col_ind: array([0, 1, 2, 3], dtype=int64)
    # X: [[1, 4, 5, 6, 5, 2, 1],
    #     [1, 4, 5, 6, 9, 24, 14, 13, 10, 6, 3, 2, 1],
    #     [1, 4, 7, 21, 19, 18, 19, 21, 7, 4, 1],
    #     [1, 4, 7, 21, 22, 27, 22, 21, 7, 4, 1]]
    # S: [11.600000000000001, 29.599999999999998, 27.2, 23.4]