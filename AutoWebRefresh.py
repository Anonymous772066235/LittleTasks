# File     :AutoWebRefresh.py
# Author   :WooChi
# Time     :2021/10/15
# Function :自动登录实验室安全中心，每隔4分钟刷新一次
# Version  :

import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from icecream import ic

import time
from selenium import webdriver  # 需pip install selenium
from selenium.webdriver.common.by import By



def refresh(sleeptime, num):
    for i in range(num):
        print(i,str((i/num)*100)+'%')
        time.sleep(sleeptime)
        driver.refresh()
    driver.close()


if __name__ == '__main__':
    # username = input('username:')
    # password = input('password:')
    username = 'WooChi'
    password = 'fWJycnm7'

    url = 'https://ca.csu.edu.cn/authserver/login?service=http%3A%2F%2F202.197.71.93%2FexamIDSLogin.php'
    driver = webdriver.Chrome()
    driver.get(url)
    driver.find_element(by=By.ID,value='username').click()  # 点击用户名输入框
    driver.find_element(by=By.ID,value='username').clear()  # 清空输入框
    driver.find_element(by=By.ID,value='username').send_keys(username)  # 自动敲入用户名

    driver.find_element(by=By.ID,value='password').click()  # 点击密码输入框
    driver.find_element(by=By.ID,value='password').clear()  # 清空输入框
    driver.find_element(by=By.ID,value='password').send_keys(password)  # 自动敲入密码

    driver.find_element(by=By.ID,value='login_submit').click()  # //*[@id="login_submit"]

    url = 'http://202.197.71.93//redir.php?catalog_id=121&object_id=2735'
    driver.get(url)
    sleeptime = 242
    num = 60
    refresh(sleeptime, num)
