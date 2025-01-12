#include <filesystem>
#include <system_error>

#include <dlfcn.h>

namespace offaxis
{
    namespace utils
    {
        std::filesystem::path abspath()
        {
            Dl_info dli;
            if (dladdr(reinterpret_cast<void *>(abspath), &dli) == 0)
                return std::filesystem::current_path() / ".";
            else
                return std::filesystem::absolute(dli.dli_fname);
        }
    }
}

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

        fp = utils::abspath().replace_filename("KBHtables80.fits");
        if (std::filesystem::exists(fp))
            return fp.string();

        throw std::system_error(std::make_error_code(std::errc::no_such_file_or_directory), fp.string());
    }
}
