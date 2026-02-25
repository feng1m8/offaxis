#include "offaxline/envs.hxx"

#include "relxill/src/LocalModel.h"

extern "C" int version_number_printed;

namespace offaxis::relxill
{
    static std::filesystem::path getFullPathTableName(const std::filesystem::path &filename)
    {
        auto env = std::getenv("OFFAXIS_TABLE_PATH");
        if (env != nullptr)
        {
            auto fp = env / filename;
            if (std::filesystem::exists(fp))
                return std::filesystem::canonical(fp);
            else
                throw std::system_error(std::make_error_code(std::errc::no_such_file_or_directory), fp.string());
        }

        env = std::getenv("XILLVER_TABLE_PATH");
        if (env != nullptr)
        {
            auto fp = env / filename;
            if (std::filesystem::exists(fp))
                return std::filesystem::canonical(fp);
            else
                throw std::system_error(std::make_error_code(std::errc::no_such_file_or_directory), fp.string());
        }

        env = std::getenv("RELXILL_TABLE_PATH");
        if (env != nullptr)
        {
            auto fp = env / filename;
            if (std::filesystem::exists(fp))
                return std::filesystem::canonical(fp);
            else
                throw std::system_error(std::make_error_code(std::errc::no_such_file_or_directory), fp.string());
        }

        auto fp = std::filesystem::current_path() / filename;
        if (std::filesystem::exists(fp))
            return std::filesystem::canonical(fp);

        fp = utils::abspath() / filename;
        if (std::filesystem::exists(fp))
            return std::filesystem::canonical(fp);

        throw std::system_error(std::make_error_code(std::errc::no_such_file_or_directory), fp.string());
    }

    static int set_xillver_table_path(T_PrimSpec prim_type)
    {
        version_number_printed = 1;

        if (prim_type == T_PrimSpec::CutoffPl)
        {
            std::string path(getFullPathTableName(XILLTABLE_FILENAME).parent_path().string());
            return putenv(utils::fotmat("RELXILL_TABLE_PATH=%s", path.c_str()).c_str());
        }

        if (prim_type == T_PrimSpec::Nthcomp)
        {
            std::string path(getFullPathTableName(XILLTABLE_NTHCOMP_FILENAME).parent_path().string());
            return putenv(utils::fotmat("RELXILL_TABLE_PATH=%s", path.c_str()).c_str());
        }

        return 0;
    }

    static xillTable *get_init_xillver_table(T_PrimSpec prim_type)
    {
        int status = set_xillver_table_path(prim_type);

        xillTable *tab = nullptr;
        get_init_xillver_table(&tab, MOD_TYPE_RELXILLLP, convertPrimSpecType(prim_type), &status);

        return tab;
    }

    int n_incl(T_PrimSpec prim_type)
    {
        if (prim_type == T_PrimSpec::CutoffPl)
        {
            static int n_incl = get_init_xillver_table(T_PrimSpec::CutoffPl)->n_incl;
            return n_incl;
        }

        if (prim_type == T_PrimSpec::Nthcomp)
        {
            static int n_incl = get_init_xillver_table(T_PrimSpec::Nthcomp)->n_incl;
            return n_incl;
        }

        return 1;
    }
}

namespace relxill
{
    [[gnu::dllexport]] void xspec_C_wrapper_eval_model(const std::vector<double> &energy, const std::vector<double> &parameter, std::vector<double> &flux, ModelName model_name)
    {
        xspec_C_wrapper_eval_model(model_name, parameter.data(), flux.data(), flux.size(), energy.data());
    }
}
