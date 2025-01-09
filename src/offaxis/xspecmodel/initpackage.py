import argparse
import errno
import os
from pathlib import Path
import shutil


def initpackage(builddir):
    builddir = Path(builddir)
    builddir.mkdir(parents=True, exist_ok=True)

    for i in builddir.iterdir():
        raise OSError(errno.ENOTEMPTY, 'Directory not empty', str(builddir))

    shutil.copyfile(Path(__file__).with_name('offaxis.dat'), builddir / 'offaxis.dat')
    shutil.copyfile(Path(__file__).with_name('liboffaxis.a'), builddir / 'liboffaxis.a')

    cwd = Path.cwd()
    os.chdir(builddir)

    os.system('initpackage offaxis offaxis.dat .')
    os.system('hmake HD_LFLAGS+="liboffaxis.a -lchealpix -fopenmp"')

    os.chdir(cwd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('builddir')
    args = parser.parse_args()

    initpackage(args.builddir)
