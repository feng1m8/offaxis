from libcpp.vector cimport vector


cdef extern from '<valarray>' namespace 'std' nogil:
    cdef cppclass valarray[T]:
        valarray() except +
        valarray(T*, size_t) except +
        valarray(size_t) except +
        T& operator[](size_t)
        size_t size()


cdef extern from 'offaxis/offaxis.h' nogil:
    void coffaxline 'offaxis::offaxline' (const valarray[double] &, const valarray[double] &, valarray[double] &) except +
    void coffaxconv 'offaxis::offaxconv' (const valarray[double] &, const valarray[double] &, valarray[double] &) except +
    void coffaxxill 'offaxis::offaxxill' (const valarray[double] &, const valarray[double] &, valarray[double] &) except +
    void coffaxxillCp 'offaxis::offaxxillCp' (const valarray[double] &, const valarray[double] &, valarray[double] &) except +


cdef xspecmodel(void (*__func__)(const valarray[double] &, const valarray[double] &, valarray[double] &) except +, str name):
    def __func__(vector[double] energy, vector[double] parameter, flux=None):
        cdef valarray[double] engs = valarray[double](energy.data(), energy.size())
        cdef valarray[double] params = valarray[double](parameter.data(), parameter.size())
        cdef valarray[double] cflux = valarray[double](energy.size() - 1)

        if flux is not None and len(flux) > 0:
            for i in range(cflux.size()):
                cflux[i] = flux[i]

        __func__(engs, params, cflux)

        if flux is None:
            return [cflux[i] for i in range(cflux.size())]

        if len(flux) > 0:
            for i in range(cflux.size()):
                flux[i] = cflux[i]
        else:
            for i in range(cflux.size()):
                flux.append(cflux[i])

    __func__ = eval('lambda energy, parameter, flux: __func__(energy, parameter, flux)', {'__func__': __func__})
    __func__.__name__ = name
    return __func__


offaxline = xspecmodel(coffaxline, 'offaxline')
offaxconv = xspecmodel(coffaxconv, 'offaxconv')
offaxxill = xspecmodel(coffaxxill, 'offaxxill')
offaxxillCp = xspecmodel(coffaxxillCp, 'offaxxillCp')
