# File     :Ransac_match.py
# Author   :WooChi
# Time     :2022/03/24
# Function :
# Version  :
import random
import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from icecream import ic

def transformation(DOPP, R, T):
    import numpy as np
    dopp_t = np.matrix(np.dot(R, DOPP.transpose()) + T)
    dopp = dopp_t.transpose()
    return np.array(dopp)


def ca_rt(DOPP, DLG,num=8,times=100,dist=1):
    import math
    import numpy as np
    repeat_ok=0
    R_ok=[]
    T_ok=[]
    for ii in range(times):
        idx = random.sample(list(np.arange(len(DLG))), num)
        ic(idx)
        Y = np.matrix(DLG[idx])
        X = np.matrix(DOPP[idx])
        mu_y = np.mean(Y, axis=0)
        mu_x = np.mean(X, axis=0)
        for i in range(len(Y)):  # 重心化
            Y[i, :] -= mu_y
        for i in range(len(X)):
            X[i, :] -= mu_x
        a = Y
        b = X
        for i in range(len(Y)):
            sin = (a[i, 0] * b[i, 1] - a[i, 1] * b[i, 0])
            cos = (a[i, 0] * b[i, 0] + a[i, 1] * b[i, 1])
        phi = math.atan(sin / cos)
        R = np.matrix([[math.cos(phi), math.sin(phi)],
                       [-math.sin(phi), math.cos(phi)]])
        T = np.matrix(mu_y.transpose() - np.dot(R, mu_x.transpose()))
        Y = np.matrix(DLG)
        X = np.matrix(DOPP)
        X2=transformation(X,R,T)
        repeat=0
        for m in range(len(X2)):
            if abs(np.around(np.linalg.norm(X2[m, :] - Y[m, :]),decimals=3))<dist:
                repeat+=1
        if repeat>repeat_ok:
            repeat_ok=repeat
            R_ok=R
            T_ok=T
        ic(repeat_ok)
    return  R_ok, T_ok,repeat_ok






def fun():
    ic()


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
