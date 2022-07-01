# File      :HistMatch.py
# Auther    :WooChi
# Time      :2022/06/29
# Version   :1.0
# Function  :

def HistMatch(img1, img2):
     """将img1直方图匹配至img2"""
     from skimage.exposure import match_histograms
     return match_histograms(img1, img2, multichannel=True)


if __name__ == '__main__':
     import cv2
     img1 = cv2.imread("D:/Project/CPP_Project/Matching/Test_Image/s1a-8bit.tif_cut_ 2.tif_cut_3.tif",0)
     img2 = cv2.imread("D:/Project/CPP_Project/Matching/Test_Image/t50-8bit.tif_cut_ 2.tif_cut_3.tif",0)
     # 变的是img1
     img3 = HistMatch(img1, img2)
     # 保存
     cv2.imwrite('D:/Project/CPP_Project/Matching/Test_Image/s1a-8bit.tif_cut_ 2.tif_cut_3.tif_HM.tif', img3)
