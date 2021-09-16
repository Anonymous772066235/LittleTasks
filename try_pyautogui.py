# File     :try_pyautogui.py
# Author   :WooChi
# Time     :2021/09/06
# Function :
# Version  :
import pyautogui as pag
pag.PAUSE=1
pag.FAILSAFE = True
from icecream import ic
ic.configureOutput(prefix='woohoo||')


if __name__=='__main__':
    x,y=pag.position()
    ic(x,y)
    ic(pag.size())
    ic(pag.onScreen(x,y))
    x=100
    y=5
    num_seconds=5
    pag.dragTo(x,y,duration=num_seconds)
    pag.click(x=x,y=y,clicks=2,interval=0.5, button='left')
    pag.scroll(5, x=x, y=y)
    pag.typewrite('hello~~~~~\n',interval=0.5)
    pag.alert('This displays some text with an OK button.')
    pag.confirm('This displays text and has an OK and Cancel button.')
    pag.prompt('This lets the user type in a string and press OK.')
    pag.screenshot('foo.png')


    import pyautogui

    print('Press Ctrl-C to quit.')
    try:
        while True:
            x, y = pyautogui.position()
            positionStr = 'X: ' + str(x).rjust(4) + ' Y: ' + str(y).rjust(4)
            print(positionStr, end='')
            print('\b' * len(positionStr), end='', flush=True)
    except KeyboardInterrupt:
        print('\n')

