import os
import errno
import shutil
import argparse
from pathlib import Path


def initpackage(builddir, force):
    builddir = Path(builddir)
    builddir.mkdir(parents=True, exist_ok=True)

    if force:
        for i in builddir.iterdir():
            if i.is_file():
                i.unlink()
            elif i.is_dir():
                shutil.rmtree(i)
    elif len(list(builddir.iterdir())) > 0:
        raise OSError(errno.ENOTEMPTY, "Directory not empty", str(builddir))

    shutil.copy(PREFIX / "offaxis.dat", builddir)

    cwd = Path.cwd()
    os.chdir(builddir)

    os.system("initpackage offaxis offaxis.dat .")
    os.system(f'hmake HD_LFLAGS+="{PREFIX}/liboffaxis_cxx.so -Wl,-rpath={PREFIX}"')

    os.chdir(cwd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("builddir")
    parser.add_argument("-f", "--force", action="store_true")
    args = parser.parse_args()

    PREFIX = Path(__file__).parent
    initpackage(args.builddir, args.force)
