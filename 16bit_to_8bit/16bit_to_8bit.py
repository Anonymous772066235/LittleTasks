# File      :16bit_to_8bit.py
# Auther    :WooChi
# Time      :2022/06/06
# Version   :2.0
# Function  :

# from YZ师兄，有改动
# 将16bit单波段影像转为8bit，在此过程中应用 “线性 2%” 拉伸

import cv2
from osgeo import gdal
import numpy as np

def transfer_16bit_to_8bit(image_path,out_path):
    #OpenCV方法，无拉伸，不考虑投影等信息
    image_16bit = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
    min_16bit = np.min(image_16bit)
    max_16bit = np.max(image_16bit)
    # image_8bit = np.array(np.rint((255.0 * (image_16bit - min_16bit)) / float(max_16bit - min_16bit)), dtype=np.uint8)
    # 或者下面一种写法
    image_8bit = np.array(np.rint(255 * ((image_16bit - min_16bit) / (max_16bit - min_16bit))), dtype=np.uint8)
    print(image_16bit.dtype)
    print('16bit dynamic range: %d - %d' % (min_16bit, max_16bit))
    print(image_8bit.dtype)
    print('8bit dynamic range: %d - %d' % (np.min(image_8bit), np.max(image_8bit)))
    cv2.imwrite(out_path,image_8bit)
    print("done")








def Liner2persent(data):
    # 线性 2% 拉伸方法
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


def U16toU8(filename,outfile):
    # GDAL方法，线性 2% 拉伸，考虑投影等信息
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

    # print(ndata.shape)
    # print(im_geotrans)
    # print(im_proj)

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
    srcpath= r'D:\Project\CPP_Project\Matching\Test_Image\s1a-iw-grd-vv-20181006t101202-20181006t101227-024014-029fab-001-pro-clip.tif'
    outpath= r'D:\Project\CPP_Project\Matching\Test_Image\s1a_8bit.tif'
    # srcpath= r'D:\Project\CPP_Project\Matching\Test_Image\t50smc_20181012t025639_b02.tif'
    # outpath= r'D:\Project\CPP_Project\Matching\Test_Image\t50-8bit.tif'
    im_data = U16toU8(srcpath,outpath)


