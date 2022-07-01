# File      :Poetry_to_password.py
# Auther    :WooChi
# Time      :2022/06/23
# Version   :1.0
# Function  :


from xpinyin import Pinyin
# from icecream import ic


def Poetry_to_password(poetry):
    # 转拼音首字母
    p = Pinyin()
    py = p.get_initials(poetry, '')
    py = py.capitalize()

    # ic(py)
    # 找到非字母字符所在位置
    symbol = []
    for i, c in enumerate(py):
        if not c.isalpha():
            symbol.append(i)

    # ic(symbol)
    # 删除非字母字符
    if symbol:
        npy = py[:symbol[0]]
        npoetry = poetry[:symbol[0]]
        for i in range(1, len(symbol)):
            npy += py[symbol[i - 1] + 1:symbol[i]].capitalize()
            npoetry += poetry[symbol[i - 1] + 1:symbol[i]]
        npy += py[symbol[- 1] + 1:].capitalize()
        npoetry += poetry[symbol[- 1] + 1:]
    else:
        npy = py
        npoetry = poetry

    # ic(npy)
    # ic(npoetry)
    # 特殊字符识别处理
    password = Special_character(npoetry, npy)

    return password


def Special_character(npoetry, npy):
    password = list(npy)
    p = Pinyin()
    py_all = p.get_pinyin(npoetry, tone_marks='marks')
    py_all = py_all.split('-')

    for i, y in enumerate(py_all):
        if y in ['líng', 'lín', 'lin', 'ling']:
            password[i] = '0'
        elif y in ['yī', 'yi']:
            password[i] = '1'
        elif y in ['èr', 'er']:
            password[i] = '2'
        elif y in ['sān', 'shān', 'san', 'shan']:
            password[i] = '3'
        elif y in ['sì', 'shì', 'si', 'shi']:
            password[i] = '4'
        elif y in ['wǔ', 'wú', 'wu', 'wǒ']:
            password[i] = '5'
        elif y in ['liù', 'niù']:
            password[i] = '6'
        elif y in ['qī', 'qi', 'xī', 'xi']:
            password[i] = '7'
        elif y in ['bā', 'ba']:
            password[i] = '8'
        elif y in ['jiǔ', 'jiu']:
            password[i] = '9'
        elif y in ['jiā']:
            password[i] = '+'
        elif y in ['jiǎn']:
            password[i] = '-'
        elif y in ['shù', 'sù']:
            password[i] = '|'

    for i, y in enumerate(npoetry):
        if y in ['云', '水', '雨', '风', '浪', '丝', '波']:
            password[i] = '~'
        elif y in ['未', '不', '惊', '叹', '无']:
            password[i] = '!'
        elif y in ['圈']:
            password[i] = '@'
        elif y in ['井', '网']:
            password[i] = '#'
        elif y in ['刀']:
            password[i] = '$'
        elif y in ['百', '白']:
            password[i] = '%'
        elif y in ['上', '高']:
            password[i] = '^'
        elif y in ['和', '兼', '且']:
            password[i] = '&'
        elif y in ['星', '乘', '日', '花', '雪']:
            password[i] = '*'
        elif y in ['下', '低', '底']:
            password[i] = '_'
        elif y in ['十']:
            password[i] = '+'
        elif y in ['等']:
            password[i] = '='
        elif y in ['大']:
            password[i] = '>'
        elif y in ['小']:
            password[i] = '<'
        elif y in ['斜', '除']:
            password[i] = '/'
        elif y in ['问', '谁', '孰', '何']:
            password[i] = '?'

    # ic(py_all)
    # ic(password)

    pw = ''
    for i, c in enumerate(password):
        pw += c
    return pw


if __name__ == '__main__':
    poetry = '晚来天欲雪，能饮一杯无，一二三四五六七八九十,和雪下网圈白'
    # poetry = '大梦谁先觉，平生我自知，加减乘除'
    # poetry = '万般皆是命，半点不由人,树斜等问'
    poetry = '不到长城非好汉，屈指行程二万。'
    # poetry = '春江潮水连海平，海上明月共潮生。' \
    #          '滟滟随波千万里，何处春江无月明！' \
    #          '江流宛转绕芳甸，月照花林皆似霰.' \
    #          '空里流霜不觉飞，汀上白沙看不见。' \
    #          '江天一色无纤尘，皎皎空中孤月轮。' \
    #          '江畔何人初见月？江月何年初照人？' \
    #          '人生代代无穷已，江月年年望相似。'\
    #          '不知江月待何人，但见长江送流水。'\
    #          '白云一片去悠悠，青枫浦上不胜愁。'\
    #          '谁家今夜扁舟子？何处相思明月楼？'\
    #          '可怜楼上月徘徊，应照离人妆镜台。'\
    #          '玉户帘中卷不去，捣衣砧上拂还来。'\
    #          '此时相望不相闻，愿逐月华流照君。'\
    #          '鸿雁长飞光不度，鱼龙潜跃水成文。'\
    #          '昨夜闲潭梦落花，可怜春半不还家。'\
    #          '江水流春去欲尽，江潭落月复西斜。'\
    #          '斜月沉沉藏海雾，碣石潇湘无限路。'\
    #          '不知乘月几人归，落月摇情满江树。'

    password = Poetry_to_password(poetry)
    print(poetry)
    print('---poetry->passward---')
    print(password)
