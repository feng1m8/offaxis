# offaxis
offaxis is a set of relativistic X-ray reflection models that allows non-axisymmetric location and velocity of X-ray source. 

- **offaxline**: Line model.
- **offaxconv**: Convolution model.
- **offaxxillCp**: Relativistic reflection model with nthcomp as incident spectrum.
- **offaxxill**: Relativistic reflection model with cutoff powerlaw as incident spectrum.

## Installation
### Required Tables
- [KBHtables80.fits](https://rec.ustc.edu.cn/share/04f322b0-0ba2-11ef-aa68-4d9b644bf13f) (required by all offaxis models)
- [xillverCp_v3.4.fits](https://rec.ustc.edu.cn/share/d2478d90-0ba2-11ef-a96a-4512f4bd416a) (required by offaxxillCp model)
- [xillver-a-Ec5.fits](https://rec.ustc.edu.cn/share/c8a256e0-0ba1-11ef-8f54-879240708725) (required by offaxxill model)

### Dependency
- xspec
- gcc or clang. (support C++17 standard)
- python>=3.8
- pkg-config
- openmp
- chealpix
- gsl
- CCfits
- cfitsio
- fftw3

### Get Source code
``` bash
git clone https://git.ustc.edu.cn/fengyuan01/offaxis.git
# or git clone https://github.com/feng1m8/offaxis.git
cd offaxis
git submodule update --init
```

### Python module
``` bash
pip install . -v
```

Try in python:
``` python
import offaxis
import xspec
xspec.Model.showList()
```

### XSPEC models
Additional step after installing Python module.
``` bash
python -m offaxis.initpackage /path/to/build/
```

Then it will be available with XSPEC.

```
lmod offaxis /path/to/build/
```

## Changelog
See [CHANGELOG.md](CHANGELOG.md).

## License
This project is distributed under [GNU GPL v3.0 or later](LICENSE).
