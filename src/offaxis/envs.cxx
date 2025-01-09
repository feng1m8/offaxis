#include <filesystem>
#include <system_error>

#include <dlfcn.h>

namespace offaxis::envs
{
    int nthreads()
    {
        auto env = std::getenv("OFFAXIS_NUM_THREADS");
        if (env == nullptr)
            return 1;
        else
            return std::atoi(env);
    }

    int nside()
    {
        auto env = std::getenv("OFFAXIS_NUM_SIDE");
        if (env == nullptr)
            return 64;
        else
            return std::atoi(env);
    }

    std::filesystem::path libpath()
    {
        Dl_info dli;
        if (dladdr((void *)libpath, &dli) == 0)
            return std::filesystem::current_path();
        else
            return std::filesystem::path(dli.dli_fname).remove_filename();
    }

    std::string kydir()
    {
        auto env = std::getenv("OFFAXIS_TABLE_PATH");
        if (env != nullptr)
        {
            auto fp = std::filesystem::path(env) / "KBHtables80.fits";
            if (std::filesystem::exists(fp))
                return fp.string();
            else
                throw std::system_error(std::make_error_code(std::errc::no_such_file_or_directory), fp.string());
        }

        env = std::getenv("KYN_TABLE_PATH");
        if (env != nullptr)
        {
            auto fp = std::filesystem::path(env) / "KBHtables80.fits";
            if (std::filesystem::exists(fp))
                return fp.string();
            else
                throw std::system_error(std::make_error_code(std::errc::no_such_file_or_directory), fp.string());
        }

        env = std::getenv("KYDIR");
        if (env != nullptr)
        {
            auto fp = std::filesystem::path(env) / "KBHtables80.fits";
            if (std::filesystem::exists(fp))
                return fp.string();
            else
                throw std::system_error(std::make_error_code(std::errc::no_such_file_or_directory), fp.string());
        }

        auto fp = std::filesystem::current_path() / "KBHtables80.fits";
        if (std::filesystem::exists(fp))
            return fp.string();

        fp = envs::libpath() / "KBHtables80.fits";
        if (std::filesystem::exists(fp))
            return fp.string();

        throw std::system_error(std::make_error_code(std::errc::no_such_file_or_directory), fp.string());
    }
}
