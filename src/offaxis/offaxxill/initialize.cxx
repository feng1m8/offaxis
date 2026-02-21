#include "offaxline/envs.hxx"

extern "C"
{
#include "relxill/src/xilltable.h"

    extern int version_number_printed;
}

namespace offaxis::relxill
{
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

        throw std::system_error(std::make_error_code(std::errc::no_such_file_or_directory), fp.string());
    }

    // static const int initialize = []()
    // {
    //     getFullPathTableName = get_full_path_table_name;
    //     version_number_printed = 1;
    //     return 0;
    // }();
}
