# File      :crop_tif_and_xml.py
# Auther    :WooChi
# Time      :2023/03/02
# Version   :1.0
# Function  :


import numpy as np
from icecream import ic
from matplotlib import pyplot as plt
from lxml import etree
import gdal
import os


def cut_img(srcpath, width_size, height_size, ColN=int(1), RowN=int(1), outpath="cut.tif"):
    """
    :param srcpath:
    :param width_size:
    :param height_size:
    :param ColN:    列向需要多少块
    :param RowN:    行向需要多少块  在这个场景下只需要一张
    :return:
    """

    in_ds = gdal.Open(srcpath)  # 读取要切的原图
    print("open tif file succeed")
    width = in_ds.RasterXSize  # 获取数据宽度
    height = in_ds.RasterYSize  # 获取数据高度
    outbandsize = in_ds.RasterCount  # 获取数据波段数
    datatype = in_ds.GetRasterBand(1).DataType
    nBands = in_ds.RasterCount  # 获取波段数

    # 读取原图中的每个波段
    in_band = []
    for n in range(nBands):
        in_band.append(in_ds.GetRasterBand(n + 1))

        # 定义切图的起始点坐标
        # 定义切图的大小（矩形框）
        col_num = int(width / width_size)  # 宽度可以分成几块
        row_num = int(height / height_size)  # 高度可以分成几块
        if (width % width_size != 0):
            col_num += 1
        if (height % height_size != 0):
            row_num += 1
        # 这边就知道我们一共是分成了多少个 如果说有多余的 那我们就让那些也自己一小块好吧
    if (ColN * RowN != 0 and ColN <= col_num and RowN <= row_num):
        col_num = ColN
        row_num = RowN

    num = 1  # 这个就用来记录一共有多少块的

    print("row_num:%d   col_num:%d" % (row_num, col_num))

    for i in range(col_num):
        for j in range(row_num):
            offset_x = i * width_size
            offset_y = j * height_size
            ## 从每个波段中切需要的矩形框内的数据(注意读取的矩形框不能超过原图大小)
            b_ysize = min(width - offset_y, height_size)
            b_xsize = min(height - offset_x, width_size)

            print("width:%d     height:%d    offset_x:%d    offset_y:%d     b_xsize:%d     b_ysize:%d" % (
            width, height, offset_x, offset_y, b_xsize, b_ysize))
            outband = []
            for n in range(nBands):
                outband.append(in_band[n].ReadAsArray(offset_y, offset_x, b_ysize, b_xsize))
            # 获取Tif的驱动，为创建切出来的图文件做准备
            gtif_driver = gdal.GetDriverByName("GTiff")

            num += 1
            # 创建切出来的要存的文件
            out_ds = gtif_driver.Create(outpath, b_ysize, b_xsize, outbandsize, datatype)
            print("create new tif file succeed")

            # 获取原图的原点坐标信息
            ori_transform = in_ds.GetGeoTransform()

            if ori_transform:
                print(ori_transform)
                print("Origin = ({}, {})".format(ori_transform[0], ori_transform[3]))
                print("Pixel Size = ({}, {})".format(ori_transform[1], ori_transform[5]))

            # 读取原图仿射变换参数值
            top_left_x = ori_transform[0]  # 左上角x坐标
            w_e_pixel_resolution = ori_transform[1]  # 东西方向像素分辨率
            top_left_y = ori_transform[3]  # 左上角y坐标
            n_s_pixel_resolution = ori_transform[5]  # 南北方向像素分辨率

            # 根据反射变换参数计算新图的原点坐标
            top_left_x = top_left_x + offset_y * w_e_pixel_resolution
            top_left_y = top_left_y + offset_x * n_s_pixel_resolution

            # 将计算后的值组装为一个元组，以方便设置
            dst_transform = (
            top_left_x, ori_transform[1], ori_transform[2], top_left_y, ori_transform[4], ori_transform[5])
            print("dst_transform = ", dst_transform)

            # 设置裁剪出来图的原点坐标
            out_ds.SetGeoTransform(dst_transform)

            # 设置SRS属性（投影信息）
            out_ds.SetProjection(in_ds.GetProjection())

            # 写入目标文件
            for n in range(nBands):
                print(outband[n])
                out_ds.GetRasterBand(n + 1).WriteArray(outband[n])
            # 将缓存写入磁盘
            out_ds.FlushCache()
            print("FlushCache succeed")
            del out_ds, outband


def zoom_bounds(bounds, width, height, width_new, height_new):
    bounds_new = np.zeros_like(bounds)
    bounds_new[0, :] = bounds[0, :]

    newRightTop = (bounds[1, :] - bounds[0, :]) * width_new / width + bounds[0, :]
    newLeftBottom = (bounds[3, :] - bounds[0, :]) * height_new / height + bounds[0, :]
    newRightBottom = (bounds[3, :] - bounds[0, :]) * height_new / height + newRightTop

    bounds_new[1, :] = newRightTop
    bounds_new[2, :] = newRightBottom
    bounds_new[3, :] = newLeftBottom

    plt.figure()
    plt.plot(bounds[:, 0], bounds[:, 1], c="r")
    plt.plot([bounds[3, 0], bounds[0, 0]], [bounds[3, 1], bounds[0, 1]], c="r")

    plt.plot(bounds_new[:, 0], bounds_new[:, 1], c="b")
    plt.plot([bounds_new[3, 0], bounds_new[0, 0]], [bounds_new[3, 1], bounds_new[0, 1]], c="b")

    plt.pause(2)

    return bounds_new


def ZY_xml_process(xmlfile, width_new, height_new, outxmlpath):
    et = etree.parse(xmlfile)

    WidthInPixels = et.find('//WidthInPixels')
    HeightInPixels = et.find('//HeightInPixels')
    width = float(WidthInPixels.text)
    height = float(HeightInPixels.text)

    if (width <= width_new or height <= height_new):
        print("Out of Bounds (width:%d height:%d)" % (width, height))
        return

    LeftTopLon = et.find('//LeftTopPoint//Longtitude')
    LeftTopLat = et.find('//LeftTopPoint//Latitude')
    RightTopLon = et.find('//RightTopPoint//Longtitude')
    RightTopLat = et.find('//RightTopPoint//Latitude')
    RightBottomLon = et.find('//RightBottomPoint//Longtitude')
    RightBottomLat = et.find('//RightBottomPoint//Latitude')
    LeftBottomLon = et.find('//LeftBottomPoint//Longtitude')
    LeftBottomLat = et.find('//LeftBottomPoint//Latitude')

    boundary = np.array([[float(LeftTopLon.text), float(LeftTopLat.text)],
                         [float(RightTopLon.text), float(RightTopLat.text)],
                         [float(RightBottomLon.text), float(RightBottomLat.text)],
                         [float(LeftBottomLon.text), float(LeftBottomLat.text)]])
    ic(boundary)

    boundary_new = zoom_bounds(boundary, width, height, width_new, height_new)

    ic(boundary_new)

    # 接下来修改xml相关属性值

    # 修改四个角点的坐标值（这部分是必须的，因为后续读取的就是这部分值）
    LeftTopLon.text = str(boundary_new[0, 0])
    LeftTopLat.text = str(boundary_new[0, 1])
    RightTopLon.text = str(boundary_new[1, 0])
    RightTopLat.text = str(boundary_new[1, 1])
    RightBottomLon.text = str(boundary_new[2, 0])
    RightBottomLat.text = str(boundary_new[2, 1])
    LeftBottomLon.text = str(boundary_new[3, 0])
    LeftBottomLat.text = str(boundary_new[3, 1])

    et.find('//RightTopPoint//Sample').text = str(width_new)
    et.find('//RightBottomPoint//Sample').text = str(width_new)
    et.find('//RightBottomPoint//Line').text = str(height_new)
    et.find('//LeftBottomPoint//Line').text = str(height_new)

    # 修改中心点以及四个角点对应的像素坐标
    center_col = int(width_new / 2)
    center_row = int(height_new / 2)
    center_lon = (boundary_new[1, 0] + boundary_new[2, 0] - boundary_new[0, 0] - boundary_new[
        3, 0]) / 2 * center_col / width_new + np.min(boundary_new[:, 0])
    center_lat = (boundary_new[0, 1] + boundary_new[1, 1] - boundary_new[2, 1] - boundary_new[
        3, 1]) / 2 * center_row / height_new + np.min(boundary_new[:, 1])

    CenterPointSample = et.find('//CenterPoint//Sample')
    CenterPointLine = et.find('//CenterPoint//Line')
    CenterPointLongtitude = et.find('//CenterPoint//Longtitude')
    CenterPointLatitude = et.find('//CenterPoint//Latitude')

    CenterPointSample.text = str(center_col)
    CenterPointLine.text = str(center_row)
    CenterPointLongtitude.text = str(center_lon)
    CenterPointLatitude.text = str(center_lat)

    # 修改宽高
    WidthInPixels.text = str(width_new)
    HeightInPixels.text = str(height_new)
    WidthInMeters = et.find('//WidthInMeters')
    HeightInMeters = et.find('//HeightInMeters')
    WidthInMeters.text = str(float(WidthInMeters.text) * width_new / width)

    HeightInMeters.text = str(float(HeightInMeters.text) * height_new / height)

    et.write(outxmlpath, encoding="utf-8")


def crop_tif_and_xml(imgPath, width_new, height_new, new_folder=""):
    (folder, imgFile) = os.path.split(imgPath)
    (name, ext) = os.path.splitext(imgFile)

    xmlFile = name + ".xml"

    xmlPath = os.path.join(folder, xmlFile)

    if (len(new_folder) == 0):
        new_folder = name[:3] + "_w" + str(width_new) + "_h" + str(height_new)

    outFolder = os.path.join(folder, new_folder)

    imgOutPath = os.path.join(outFolder, imgFile)
    xmlOutPath = os.path.join(outFolder, xmlFile)

    ic(outFolder)
    ic(imgOutPath)
    ic(xmlOutPath)

    if not os.path.exists(outFolder):
        os.makedirs(outFolder)  # makedirs 创建文件时如果路径不存在会创建这个路径

    if (os.path.exists(imgPath)):
        print(ext, " processing...")
        cut_img(imgPath, int(width_new), int(height_new), outpath=imgOutPath)
    else:
        print("no img file")
    if (os.path.exists(xmlPath) and (name[:2] == "zy" or name[:2] == "ZY")):
        print("xml processing...")
        ZY_xml_process(xmlPath, width_new, height_new, xmlOutPath)
    else:
        print("no xml file")

    print("done")


if __name__ == '__main__':
    imgPath = r"T:\ProjectData\ZY302_2YW\ZY302_2YW\zy302a_bwd_016497_003126_20190520112310_01_sec_0001_1905228482.tif"
    width_new = 3000
    height_new = 3000

    crop_tif_and_xml(imgPath, width_new, height_new)
