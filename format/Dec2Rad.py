# File     :Dec2Rad.py
# Author   :WooChi
# Time     :2021/11/24
# Function :经纬度转换
# Version  :


import time
from icecream import ic

# [十进制小数]
# 转换为[度分秒]
def LatLng_Dec2Rad(decNum):
    NumIntegral = int(decNum)  # 整数部分
    NumDecimal = decNum - NumIntegral  # 小数部分

    tmp = NumDecimal * 3600
    degree = NumIntegral  # 度
    minute = int(tmp // 60)  # 分
    second = tmp - minute * 60  # 秒 tmp%3600

    return degree, minute, second


# [度分秒] 转换为 [十进制小数]
def LatLng_Rad2Dec(d, m, s):
    decNum = d + m / 60.0 + s / 3600.0
    return decNum


if __name__ == "__main__":
    lng_decNum =-113.211  # 要转换的经度

    t = LatLng_Dec2Rad(lng_decNum)
    print(t)

    k = LatLng_Rad2Dec(t[0], t[1], t[2])
    print(k)

