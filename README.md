# offaxis
offaxis is a set of relativistic X-ray reflection models that allows non-axisymmetric location and velocity of X-ray source. 

- **offaxline**: Line model.
- **offaxconv**: Convolution model.
- **offaxxillCp**: Relativistic reflection model with nthcomp as incident spectrum.
- **offaxxill**: Relativistic reflection model with cutoff powerlaw as incident spectrum.

<!-- ## Installation -->
## Required Tables
- [KBHtables80.fits](https://rec.ustc.edu.cn/share/04f322b0-0ba2-11ef-aa68-4d9b644bf13f) (required by all offaxis models)
- [xillverCp_v3.4.fits](https://rec.ustc.edu.cn/share/d2478d90-0ba2-11ef-a96a-4512f4bd416a) (required by offaxxillCp model)
- [xillver-a-Ec5.fits](https://rec.ustc.edu.cn/share/c8a256e0-0ba1-11ef-8f54-879240708725) (required by offaxxill model)

<!-- ### Dependency
- Python>=3.7 (Recommended Python>=3.10)
- C and C++ compiler GCC or Clang.
- pkg-config
- openmp
- CCfits
- chealpix
- gsl
- fftw3 -->

<!-- ### Installing -->
<!-- 1. First, clone this project to local computer:

       git clone --recurse-submodules https://github.com/feng1m8/offaxline.git -->

<!-- 2. It is recommended to install with pip:

       cd offaxline
       pip install . -v

    Now it is available as PyXspec models.

       bash:~$ python
       >>> import offaxline
       >>> import xspec
       >>> xspec.Model.showList()
	    Additive Models:
	   ...
	   offaxline*  offaxxill*  offaxxillCp*
	   ...

		Convolution Models:
	   ...
	   offaxconv*
	   ...

3. Additional step required to use it in XSPEC.

       python -m offaxline.initpackage /path/to/build/ --force

    Then it will be available with XSPEC.

       bash:~$ xspec
       XSPEC> load /path/to/build/liboffaxline.so -->

## Changelog
See [CHANGELOG.md](CHANGELOG.md).

## License
This project is distributed under [GNU GPL v3.0 or later](LICENSE).

<!-- ## Reference

    @ARTICLE{offaxline,
        author = {{Feng}, Yuan and {Yuan}, Ye-Fei and {Zhang}, Shuang-Nan},
        title = "{Moving Corona and the Line Profile of the Relativistic Broad Iron Emission Line}",
        journal = {\apj},
        year = 2023,
        month = sep,
        volume = {955},
        number = {1},
        pages = {53},
        doi = {10.3847/1538-4357/acedff},
    } -->
