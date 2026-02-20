#include <cstdarg>
#include <filesystem>

#include <dlfcn.h>

namespace offaxis
{
    namespace utils
    {
        std::filesystem::path abspath()
        {
            Dl_info dli;
            if (dladdr(reinterpret_cast<void *>(abspath), &dli) == 0)
                return std::filesystem::current_path();
            else
                return std::filesystem::absolute(dli.dli_fname).parent_path();
        }

        std::string fotmat(const char *format, ...)
        {
            std::va_list args;
            va_start(args, format);
            std::string buffer(std::vsnprintf(nullptr, 0, format, args), '\0');
            va_end(args);

            va_start(args, format);
            std::vsnprintf(buffer.data(), buffer.size() + 1, format, args);
            va_end(args);

            return buffer;
        }
    }

    namespace envs
    {
        int nthreads()
        {
            auto env = std::getenv("OFFAXIS_NUM_THREADS");
            if (env == nullptr)
                return 1;
            else
                return std::atoi(env);
        }
    }
}
