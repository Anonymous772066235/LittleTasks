# File     :Unsplash.py
# Author   :WooChi
# Time     :2021/09/06
# Function :
# Version  :
# import requests
# import os
# from bs4 import BeautifulSoup
# import re
# import bs4
# import json
# url = "https://unsplash.com/napi/photos?page="
# base = "&per_page=12"
# def getHTML(url):
#     r = requests.get(url)
#     jd = json.loads(r.text)
#     list = []
#     for i in jd:
#         t = i['urls']["raw"]
#         list.append(t)
#     root = "F:\\Chrome Downloads\\Pic"
#     global count
#     for i in list:
#         print(i)
#         path = root + str(count) + '.jpg'
#         count = count + 1
#         if not os.path.exists(root):
#             os.mkdir(root)
#         if not os.path.exists(path):
#             r = requests.get(i)
#
#             with open(path, 'wb') as f:
#                 f.write(r.content)
# count = 0
# for i in range(1,3):
#     u = url + str(i)+base
#     print(u)
#     getHTML(u)

# -*- coding:UTF-8 -*-
import requests, json, time
from contextlib import closing
from icecream import ic
ic.configureOutput(prefix='woohoo||')
class get_photos(object):
	def __init__(self):
		self.photos_id = []
		self.download_server = 'https://unsplash.com/photos/xxx/download?force=true'
		self.target = 'https://unsplash.com/napi/photos?page=xxx&per_page=12'

	'''
	函数说明：获取图片ID
	Parameters:
		page --页数
	Returns:
		None
	'''

	def get_ids(self, page):
		target = self.target.replace('xxx', str(page))
		req = requests.get(url = target)
		html = json.loads(req.text)
		for each in html:
			self.photos_id.append(each['id'])
		time.sleep(1)
		for i in range(5):	#获取6页图片的id
			page = page + 1
			next_page = self.target.replace('xxx', str(page))
			req = requests.get(url = next_page)
			html = json.loads(req.text)
			for each in html:
				self.photos_id.append(each['id'])
			time.sleep(1)

	'''
	函数说明：图片下载
	Parameters:
		photo_id --图片id
		filename --图片存储名
	Returns:
		None
	'''
	def download(self, photo_id, filename):
		ic(filename)
		target = self.download_server.replace('xxx', photo_id)
		ic(target)
		with closing(requests.get(url = target, stream = True)) as r:
			with open('%s.jpg' % filename, 'ab+') as f:
				for chunk in r.iter_content(chunk_size = 1024):
					if chunk:
						f.write(chunk)
						f.flush()

if __name__ == '__main__':
	gp = get_photos()
	print('获取图片链接中：')
	gp.get_ids(1)
	print('图片下载中。。。')
	for i in range(len(gp.photos_id)):
		print('正在下载第%d张图片' % (i+1))
		gp.download(gp.photos_id[i], 'F:\\Chrome Downloads\\Pic\\'+gp.photos_id[i])	#使用图片id作为图片的存储名
