import os
import errno
import shutil
import argparse
import subprocess
from pathlib import Path


def initpackage(builddir, force):
    builddir = Path(builddir)
    builddir.mkdir(parents=True, exist_ok=True)

    if force:
        for i in builddir.iterdir():
            if i.is_dir(follow_symlinks=False):
                shutil.rmtree(i)
            else:
                i.unlink()
    elif not list(builddir.iterdir()):
        raise OSError(errno.ENOTEMPTY, "Directory not empty", str(builddir))

    PREFIX = Path(__file__).parent
    shutil.copy2(PREFIX / "offaxis.dat", builddir)

    cwd = Path.cwd()
    os.chdir(builddir)

    subprocess.run(["initpackage", "offaxis", "offaxis.dat", "."])
    subprocess.run(["hmake", f'HD_LFLAGS+="{PREFIX}/liboffaxis_cxx.so"', f"-Wl,-rpath={PREFIX}"])

    os.chdir(cwd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("builddir")
    parser.add_argument("-f", "--force", action="store_true")
    args = parser.parse_args()

    initpackage(args.builddir, args.force)
