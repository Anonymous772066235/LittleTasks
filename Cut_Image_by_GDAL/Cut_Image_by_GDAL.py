# File      :Cut_Image_by_GDAL.py
# Auther    :WooChi
# Time      :2022/06/12
# Version   :1.0
# Function  :


import gdal


def Cut(srcpath,width_size,height_size,offset_x = 0,offset_y = 0):
    """
    :param srcpath:
    :param width_size:
    :param height_size:
    :param offset_x: 起始坐标，默认为0
    :param offset_y:
    :return: None
    """
    in_ds = gdal.Open(srcpath)              # 读取要切的原图
    print("open tif file succeed")
    width = in_ds.RasterXSize                         # 获取数据宽度
    height = in_ds.RasterYSize                        # 获取数据高度
    outbandsize = in_ds.RasterCount                   # 获取数据波段数
    datatype = in_ds.GetRasterBand(1).DataType
    nBands=in_ds.RasterCount                       # 获取波段数

    # 读取原图中的每个波段
    in_band=[]
    for n in range(nBands):
        in_band.append(in_ds.GetRasterBand(n+1))

    # 定义切图的起始点坐标
    # 定义切图的大小（矩形框）
    col_num = int(width / width_size)  #宽度可以分成几块
    row_num = int(height / height_size) #高度可以分成几块
    if(width % width_size!= 0):
        col_num += 1
    if(height % height_size != 0):
        row_num += 1
    #这边就知道我们一共是分成了多少个 如果说有多余的 那我们就让那些也自己一小块好吧

    num = 1    #这个就用来记录一共有多少块的

    print("row_num:%d   col_num:%d" %(row_num,col_num))

    for i in range(col_num):
        for j in range(row_num):
            offset_x = i * width_size
            offset_y = j * height_size
            ## 从每个波段中切需要的矩形框内的数据(注意读取的矩形框不能超过原图大小)
            b_ysize = min(width - offset_y, height_size)
            b_xsize = min(height - offset_x, width_size)

            print("width:%d     height:%d    offset_x:%d    offset_y:%d     b_xsize:%d     b_ysize:%d" %(width,height,offset_x,offset_y, b_xsize, b_ysize))
            # print("\n")
            outband=[]
            for n in range(nBands):
                outband.append(in_band[n].ReadAsArray(offset_y, offset_x, b_ysize, b_xsize))
            # 获取Tif的驱动，为创建切出来的图文件做准备
            gtif_driver = gdal.GetDriverByName("GTiff")
            file = srcpath+"_cut_%1d" % num+".tif"
            num += 1
            # 创建切出来的要存的文件
            out_ds = gtif_driver.Create(file, b_ysize, b_xsize, outbandsize, datatype)
            print("create new tif file succeed")

            # 获取原图的原点坐标信息
            ori_transform = in_ds.GetGeoTransform()

            if ori_transform:
                    print (ori_transform)
                    print("Origin = ({}, {})".format(ori_transform[0], ori_transform[3]))
                    print("Pixel Size = ({}, {})".format(ori_transform[1], ori_transform[5]))

            # 读取原图仿射变换参数值
            top_left_x = ori_transform[0]  # 左上角x坐标
            w_e_pixel_resolution = ori_transform[1] # 东西方向像素分辨率
            top_left_y = ori_transform[3] # 左上角y坐标
            n_s_pixel_resolution = ori_transform[5] # 南北方向像素分辨率

            # 根据反射变换参数计算新图的原点坐标
            top_left_x = top_left_x + offset_y * w_e_pixel_resolution
            top_left_y = top_left_y + offset_x * n_s_pixel_resolution

            # 将计算后的值组装为一个元组，以方便设置
            dst_transform = (top_left_x, ori_transform[1], ori_transform[2], top_left_y, ori_transform[4], ori_transform[5])
            print("dst_transform = ",dst_transform)

            # 设置裁剪出来图的原点坐标
            out_ds.SetGeoTransform(dst_transform)

            # 设置SRS属性（投影信息）
            out_ds.SetProjection(in_ds.GetProjection())

            # 写入目标文件
            for n in range(nBands):
                print(outband[n])
                out_ds.GetRasterBand(n+1).WriteArray(outband[n])
            # 将缓存写入磁盘
            out_ds.FlushCache()
            print("FlushCache succeed")
            del out_ds,outband


if __name__ == '__main__':
    # srcpath=r'D:\Project\CPP_Project\Matching\Test_Image\s1a-8bit.tif'
    # srcpath=r'D:\Project\CPP_Project\Matching\Test_Image\t50-8bit.tif'
    # srcpath=r'D:\Project\CPP_Project\Matching\Test_Image\s1a-8bit.tif_cut_ 2.tif'
    # srcpath=r'D:\Project\CPP_Project\Matching\Test_Image\t50-8bit.tif_cut_ 2.tif'
    # srcpath=r'D:\Project\CPP_Project\Matching\Test_Image\s1a-8bit.tif_cut_ 2.tif_cut_3.tif'
    srcpath=r'D:\Project\CPP_Project\Matching\Test_Image\t50-8bit.tif_cut_ 2.tif_cut_3.tif'
    Cut(srcpath,int(1000),int(1000) )
    pass
