# File     :creat_new_py.py
# Author   :WooChi
# Time     :2021/09/08
# Function :
# Version  :

import pyautogui as pag
pag.PAUSE=0.5
pag.FAILSAFE = True
from icecream import ic
ic.configureOutput(prefix='woohoo||')

def creat_new_py(py_name):


    pag.click(x=55,y=12,clicks=1,button='left')
    pag.hotkey('Alt','Insert')



    pag.hotkey('Down')
    pag.hotkey('Down')
    pag.hotkey('Down')
    pag.hotkey('Down')
    pag.hotkey('Enter')
    pag.typewrite(py_name)
    pag.hotkey('Enter','Enter')



if __name__=='__main__':
    py_name='LSDtest1104'
    creat_new_py(py_name)

    # while True:
    #     print(pag.position())
