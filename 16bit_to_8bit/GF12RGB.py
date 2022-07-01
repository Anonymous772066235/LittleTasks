# from YZ师兄，有改动
# 将16bit单波段影像转为8bit，在此过程中应用 “线性 2%” 拉伸

import cv2
from osgeo import gdal
import numpy as np
from icecream import ic

def Liner2persent(data):
    mdata = data.copy()
    rows, cols = data.shape[0:2]
    counts = rows*cols

    mdata = mdata.reshape(counts, 1)
    tdata = np.sort(mdata, axis=0)
    cutmin = tdata[int(counts*0.01), 0]
    cutmax = tdata[int(counts*0.99), 0]

    ndata = 255.0*(data.astype(np.float32) - cutmin)/float(cutmax-cutmin)
    ndata[data < cutmin] = 0
    ndata[data > cutmax] = 255
    return ndata

# def U16_4Channels2U8_3channels(filename):
#     dataset = gdal.Open(filename)
#     im_width = dataset.RasterXSize  #栅格矩阵的列数
#     im_height = dataset.RasterYSize  #栅格矩阵的行数
#     im_geotrans = dataset.GetGeoTransform() #仿射矩阵
#     im_proj = dataset.GetProjection() #地图投影信息
#     im_data = dataset.ReadAsArray(0,0,im_width,im_height) #将数据写成数组，对应栅格矩阵
#     r = (Liner2persent(im_data[0,:,:])).astype(np.uint8)
#     g = (Liner2persent(im_data[1,:,:])).astype(np.uint8)
#     b = (Liner2persent(im_data[2,:,:])).astype(np.uint8)
#     return np.stack((r, g, b), axis=2)


def U16toU8(filename,outfile):
    dataset = gdal.Open(filename)
    im_width = dataset.RasterXSize  #栅格矩阵的列数
    im_height = dataset.RasterYSize  #栅格矩阵的行数
    im_geotrans = dataset.GetGeoTransform() #仿射矩阵
    im_proj = dataset.GetProjection() #地图投影信息
    nBands=dataset.RasterCount

    if nBands!=1:
        print("仅针对单波段影像，而该影响为%d波段影像"%nBands)
        return -1

    im_data = dataset.ReadAsArray(0,0,im_width,im_height) #将数据写成数组，对应栅格矩阵

    ndata = Liner2persent(im_data).astype(np.uint8)
    print(ndata)



    print(ndata.shape)
    print(im_geotrans)
    print(im_proj)

    # 保存影像
    gtif_driver = gdal.GetDriverByName("GTiff")
    # 创建切出来的要存的文件
    out_ds = gtif_driver.Create(outfile, im_width, im_height, nBands  ,gdal.GDT_Byte )

    print("create new tif file succeed")

    out_ds.SetGeoTransform(im_geotrans)
    out_ds.SetProjection(im_proj)

    # 写入目标文件
    out_ds.GetRasterBand(1).WriteArray(ndata)
    # 将缓存写入磁盘
    out_ds.FlushCache()
    print("FlushCache succeed")
    del out_ds, im_data

    return ndata

if __name__ == "__main__":
    srcpath= r'D:\Project\CPP_Project\Matching\Test_Image\ReferenceImage\s1a-iw-grd-vv-20181006t101202-20181006t101227-024014-029fab-001-pro-clip.tif'
    outpath = r'D:\Project\CPP_Project\Matching\Test_Image\ReferenceImage\s1a_8bit.tif'
    # im_data = U16_4Channels2U8_3channels(filename)
    im_data = U16toU8(srcpath,outpath)
    cv2.namedWindow('test', 0)
    cv2.imwrite('D:/Project/CPP_Project/Matching/Test_Image/InputImage/sla_py.tif',im_data)
    cv2.imshow('test', im_data)
    cv2.waitKey()
