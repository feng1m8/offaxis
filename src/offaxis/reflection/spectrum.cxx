#include <functional>
#include <memory>

#include <stdlib.h>

#include "offaxis/convolve/convolve.hxx"
#include "offaxis/envs.hxx"
#include "offaxis/parameter.hxx"
#include "spectrum.hxx"

extern "C" int version_number_printed;

namespace offaxis::relxill
{
    Spectrum::Spectrum(const std::valarray<double> &parameter, T_PrimSpec prim_type) : ect(this->param.ect)
    {
        using namespace parameter::offaxxillCp;

        this->param.gam = parameter[gamma];
        this->param.afe = parameter[Afe];
        this->param.lxi = parameter[logxi];
        this->param.ect = parameter[Ecut];
        this->param.incl = parameter[Incl];
        this->param.dens = parameter[logN];
        this->param.frac_pl_bb = 0.0;
        this->param.kTbb = 0.0;
        this->param.prim_type = convertPrimSpecType(prim_type);
        this->param.model_type = MOD_TYPE_RELXILLLP;
    }

    double Spectrum::norm(double ener_shift_observer_source) const
    {
        return calc_xillver_normalization_change(1.0 / ener_shift_observer_source, &this->param);
    }

    std::valarray<double> Spectrum::angdep(const std::valarray<double> &energy, const std::valarray<double> &hist, const std::valarray<double> &dist) const
    {
        int status = EXIT_SUCCESS;

        std::unique_ptr<xillSpec, std::function<void(xillSpec *)>> xill_spec(get_xillver_spectra_table(&this->param, &status), free_xill_spec);
        std::valarray<double> xill_flux(xill_spec->n_ener);
        calc_xillver_angdep(std::begin(xill_flux), xill_spec.get(), std::begin(dist), &status);

        relxill::convolveSpectrumFFTNormalized(energy, hist, xill_spec->ener, xill_flux);

        return xill_flux;
    }

    std::valarray<double> Spectrum::nthcomp(const std::valarray<double> &ener, double ener_shift_source_obs) const
    {
        static const std::valarray<double> energy_norm(utils::geomspace(EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION, N_ENER_COARSE + 1));

        int status = EXIT_SUCCESS;

        double prim_spec_source[N_ENER_COARSE];
        calc_primary_spectrum(prim_spec_source, std::begin(energy_norm), N_ENER_COARSE, &this->param, &status);
        double norm = 2.0 / calcNormWrtXillverTableSpec(prim_spec_source, std::begin(energy_norm), N_ENER_COARSE, nullptr);

        std::valarray<double> flux(ener.size() - 1);
        calc_primary_spectrum(std::begin(flux), std::begin(ener), flux.size(), &this->param, &status, ener_shift_source_obs);

        return norm * flux;
    }

    int n_incl(T_PrimSpec prim_type)
    {
        xillTable *tab = nullptr;
        int status = EXIT_SUCCESS;

        get_init_xillver_table(&tab, MOD_TYPE_RELXILLLP, convertPrimSpecType(prim_type), &status);

        return tab->n_incl;
    }

    char *get_full_path_table_name(const char *filename, int *status)
    {
        auto env = std::getenv("OFFAXIS_TABLE_PATH");
        if (env != nullptr)
        {
            auto fp = std::filesystem::path(env) / filename;
            if (std::filesystem::exists(fp))
            {
                char *fullfilename = new char[1 + fp.string().size()]{'\0'};
                fp.string().copy(fullfilename, fp.string().size());
                return fullfilename;
            }
            else
                throw std::system_error(std::make_error_code(std::errc::no_such_file_or_directory), fp.string());
        }

        env = std::getenv("RELXILL_TABLE_PATH");
        if (env != nullptr)
        {
            auto fp = std::filesystem::path(env) / filename;
            if (std::filesystem::exists(fp))
            {
                char *fullfilename = new char[1 + fp.string().size()]{'\0'};
                fp.string().copy(fullfilename, fp.string().size());
                return fullfilename;
            }
            else
                throw std::system_error(std::make_error_code(std::errc::no_such_file_or_directory), fp.string());
        }

        auto fp = std::filesystem::current_path() / filename;
        if (std::filesystem::exists(fp))
        {
            char *fullfilename = new char[1 + fp.string().size()]{'\0'};
            fp.string().copy(fullfilename, fp.string().size());
            return fullfilename;
        }

        fp = utils::abspath().replace_filename(filename);
        if (std::filesystem::exists(fp))
        {
            char *fullfilename = new char[1 + fp.string().size()]{'\0'};
            fp.string().copy(fullfilename, fp.string().size());
            return fullfilename;
        }
        else
            throw std::system_error(std::make_error_code(std::errc::no_such_file_or_directory), fp.string());

        return nullptr;
    }

    static const int init = []()
    {
        getFullPathTableName = get_full_path_table_name;
        version_number_printed = 1;
        return 0;
    }();
}
