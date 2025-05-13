from libcpp.vector cimport vector


cdef extern from '<valarray>' namespace 'std' nogil:
    cdef cppclass valarray[T]:
        valarray() except +
        valarray(T*, size_t) except +
        valarray(size_t) except +
        T& operator[](size_t)
        size_t size()


cdef extern from 'offaxis/offaxis.h':
    void C_offaxline 'offaxis::offaxline' (const valarray[double] &, const valarray[double] &, valarray[double] &) except +
    void C_offaxconv 'offaxis::offaxconv' (const valarray[double] &, const valarray[double] &, valarray[double] &) except +
    void C_offaxxill 'offaxis::offaxxill' (const valarray[double] &, const valarray[double] &, valarray[double] &) except +
    void C_offaxxillCp 'offaxis::offaxxillCp' (const valarray[double] &, const valarray[double] &, valarray[double] &) except +


cdef xspecmodel(void (*func)(const valarray[double] &, const valarray[double] &, valarray[double] &) except +):
    def model(vector[double] energy, vector[double] parameter, flux=None):
        cdef valarray[double] engs = valarray[double](energy.data(), energy.size())
        cdef valarray[double] param = valarray[double](parameter.data(), parameter.size())
        cdef valarray[double] cflux = valarray[double](energy.size() - 1)

        if flux is not None and len(flux) > 0:
            for i in range(cflux.size()):
                cflux[i] = flux[i]

        func(engs, param, cflux)

        if flux is None:
            return [cflux[i] for i in range(cflux.size())]

        if len(flux) > 0:
            for i in range(cflux.size()):
                flux[i] = cflux[i]
        else:
            for i in range(cflux.size()):
                flux.append(cflux[i])

    return model


offaxline = eval('lambda energy, parameter, flux: C_offaxline(energy, parameter, flux)', {'C_offaxline': xspecmodel(&C_offaxline)})
offaxconv = eval('lambda energy, parameter, flux: C_offaxconv(energy, parameter, flux)', {'C_offaxconv': xspecmodel(&C_offaxconv)})
offaxxill = eval('lambda energy, parameter, flux: C_offaxxill(energy, parameter, flux)', {'C_offaxxill': xspecmodel(&C_offaxxill)})
offaxxillCp = eval('lambda energy, parameter, flux: C_offaxxillCp(energy, parameter, flux)', {'C_offaxxillCp': xspecmodel(&C_offaxxillCp)})

offaxline.__name__ = 'offaxline'
offaxconv.__name__ = 'offaxconv'
offaxxill.__name__ = 'offaxxill'
offaxxillCp.__name__ = 'offaxxillCp'
