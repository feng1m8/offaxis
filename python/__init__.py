from .offaxis import *


def lmod():
    from xspec import AllModels
    from pathlib import Path

    info = Path(__file__).with_name('offaxis.dat').read_text().split('\n\n')
    info = [tuple(i.splitlines()[1:]) for i in info]

    AllModels.addPyMod(offaxline, info[0], 'add')
    AllModels.addPyMod(offaxconv, info[1], 'con')
    AllModels.addPyMod(offaxxill, info[2], 'add')
    AllModels.addPyMod(offaxxillCp, info[3], 'add')


try:
    lmod()
except:
    pass


url = {
    'KBHtables': 'https://rec.ustc.edu.cn/share/04f322b0-0ba2-11ef-aa68-4d9b644bf13f',
    'xillver': 'https://rec.ustc.edu.cn/share/c8a256e0-0ba1-11ef-8f54-879240708725',
    'xillverCp': 'https://rec.ustc.edu.cn/share/d2478d90-0ba2-11ef-a96a-4512f4bd416a',
}
