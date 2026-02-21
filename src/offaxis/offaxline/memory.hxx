#ifndef OFFAXLINE_MEMORY_HXX
#define OFFAXLINE_MEMORY_HXX

#include <algorithm>
#include <deque>
#include <functional>
#include <string>

namespace offaxis
{
    namespace envs
    {
        inline std::size_t cache_size()
        {
            auto env = std::getenv("OFFAXIS_CACHE_SIZE");
            if (env == nullptr)
                return 0;
            else
                return std::stoull(env);
        }
    }

    template <typename R, typename... Args>
    class Memory
    {
    private:
        const std::function<R(Args...)> func;
        std::deque<std::pair<std::tuple<typename std::remove_reference<Args>::type...>, R>> data;

    public:
        std::size_t max_size;

        Memory(R (*f)(Args...)) : func(f), max_size(0) {};

        R operator()(Args... args)
        {
            if (this->max_size == 0)
            {
                return func(args...);
            }

            auto params = std::tie(args...);

            auto i = std::find_if(
                this->data.cbegin(), this->data.cend(),
                [&](auto &i)
                {
                    return i.first == params;
                });

            if (i == this->data.cend())
            {
                auto res = std::make_pair(params, func(args...));

                if (this->data.size() >= this->max_size)
                    this->data.resize(this->max_size - 1);

                return this->data.emplace_front(res).second;
            }

            return i->second;
        }
    };
}

#endif
