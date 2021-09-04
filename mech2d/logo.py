from mech2d.utils import box_center
def logo():

    logo_list=     ['                     _     ___     _ ',
                    '                    | |   |__ \   | |',
                    ' _ __ ___   ___  ___| |__    ) |__| |',
                    "| '_ ` _ \ / _ \/ __| '_ \  / // _` |",
                    '| | | | | |  __/ (__| | | |/ /| (_| |',
                    '|_| |_| |_|\___|\___|_| |_|____\__,_|']

    box_center(ch='_',fill='_',sp=' ')
    box_center(ch=' ',fill=' ',sp='|')
    for lstr in logo_list:
        box_center(ch=lstr,fill=' ',sp='|')
    box_center(ch=' ',fill=' ',sp='|')
    box_center(ch='version 1.0.0',fill=' ',sp='|')
    box_center(ch='_',fill='_',sp='|')
