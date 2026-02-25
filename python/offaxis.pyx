from libcpp.vector cimport vector


cdef extern from "<valarray>" namespace "std" nogil:
    cdef cppclass valarray[T]:
        valarray() except +
        valarray(T*, size_t) except +
        valarray(size_t) except +
        T& operator[](size_t)
        size_t size()


cdef extern from "offaxis/offaxis.h" nogil:
    void coffaxline "offaxis::offaxline" (const valarray[double] &, const valarray[double] &, valarray[double] &) except +
    void coffaxconv "offaxis::offaxconv" (const valarray[double] &, const valarray[double] &, valarray[double] &) except +
    void coffaxxill "offaxis::offaxxill" (const valarray[double] &, const valarray[double] &, valarray[double] &) except +
    void coffaxxillCp "offaxis::offaxxillCp" (const valarray[double] &, const valarray[double] &, valarray[double] &) except +


cdef xspecmodel(void (*F)(const valarray[double] &, const valarray[double] &, valarray[double] &) except +, str name):
    def __func__(vector[double] energy, vector[double] parameter, flux=None):
        cdef valarray[double] engs = valarray[double](energy.data(), energy.size())
        cdef valarray[double] params = valarray[double](parameter.data(), parameter.size())
        cdef valarray[double] cflux = valarray[double](energy.size() - 1)

        if flux is not None and len(flux) > 0:
            for i in range(cflux.size()):
                cflux[i] = flux[i]

        F(engs, params, cflux)

        if flux is None:
            return [cflux[i] for i in range(cflux.size())]

        if len(flux) > 0:
            for i in range(cflux.size()):
                flux[i] = cflux[i]
        else:
            for i in range(cflux.size()):
                flux.append(cflux[i])

    __func__ = eval("lambda energy, parameter, flux=None: __func__(energy, parameter, flux)", {"__func__": __func__})
    __func__.__name__ = name
    return __func__


offaxline = xspecmodel(coffaxline, "offaxline")
offaxconv = xspecmodel(coffaxconv, "offaxconv")
offaxxill = xspecmodel(coffaxxill, "offaxxill")
offaxxillCp = xspecmodel(coffaxxillCp, "offaxxillCp")


cdef extern from "relxill/src/ModelInfo.h" nogil:
    """
namespace relxill
{
    void xspec_C_wrapper_eval_model(const std::vector<double> &, const std::vector<double> &, std::vector<double> &, ModelName);
}
    """

    cdef enum class ModelName:
        relline
        relline_lp
        relline_ext
        relconv
        relconv_lp
        relxill
        relxillNS
        relxillCO
        relxillCp
        relxilllp
        relxill_ext_ecut
        relxilllpCp
        relxill_ext
        relxilllpion
        relxilllpionCp
        relxilllpAlpha
        xillver
        xillverCp
        xillverNS
        xillverCO
        xillverD
        relxillD
        relxilllpD
        relxillBB
        relxill_jedsad

    void c_xspec_C_wrapper_eval_model "relxill::xspec_C_wrapper_eval_model"(const vector[double] &, const vector[double] &, vector[double] &, ModelName) except +


cdef xspec_C_wrapper_eval_model(ModelName model_name):
    def __func__(vector[double] energy, vector[double] parameter, flux=None):
        cdef vector[double] cflux = vector[double](energy.size() - 1)

        if flux is not None and len(flux) > 0:
            for i in range(cflux.size()):
                cflux[i] = flux[i]

        c_xspec_C_wrapper_eval_model(energy, parameter, cflux, model_name)

        if flux is None:
            return [cflux[i] for i in range(cflux.size())]

        if len(flux) > 0:
            for i in range(cflux.size()):
                flux[i] = cflux[i]
        else:
            for i in range(cflux.size()):
                flux.append(cflux[i])

    return __func__


relline = xspec_C_wrapper_eval_model(ModelName.relline)
relline_lp = xspec_C_wrapper_eval_model(ModelName.relline_lp)
relline_ext = xspec_C_wrapper_eval_model(ModelName.relline_ext)
relconv = xspec_C_wrapper_eval_model(ModelName.relconv)
relconv_lp = xspec_C_wrapper_eval_model(ModelName.relconv_lp)
relxill = xspec_C_wrapper_eval_model(ModelName.relxill)
relxillNS = xspec_C_wrapper_eval_model(ModelName.relxillNS)
relxillCO = xspec_C_wrapper_eval_model(ModelName.relxillCO)
relxillCp = xspec_C_wrapper_eval_model(ModelName.relxillCp)
relxilllp = xspec_C_wrapper_eval_model(ModelName.relxilllp)
relxill_ext_ecut = xspec_C_wrapper_eval_model(ModelName.relxill_ext_ecut)
relxilllpCp = xspec_C_wrapper_eval_model(ModelName.relxilllpCp)
relxill_ext = xspec_C_wrapper_eval_model(ModelName.relxill_ext)
relxilllpAlpha = xspec_C_wrapper_eval_model(ModelName.relxilllpAlpha)
xillver = xspec_C_wrapper_eval_model(ModelName.xillver)
xillverCp = xspec_C_wrapper_eval_model(ModelName.xillverCp)
xillverNS = xspec_C_wrapper_eval_model(ModelName.xillverNS)
xillverCO = xspec_C_wrapper_eval_model(ModelName.xillverCO)
relxillBB = xspec_C_wrapper_eval_model(ModelName.relxillBB)
relxill_jedsad = xspec_C_wrapper_eval_model(ModelName.relxill_jedsad)
