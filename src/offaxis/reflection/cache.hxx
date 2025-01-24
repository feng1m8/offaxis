#include <algorithm>
#include <deque>
#include <functional>

namespace offaxis
{
    namespace envs
    {
        inline int cache_size()
        {
            auto env = std::getenv("OFFAXIS_CACHE_SIZE");
            if (env == nullptr)
                return 1;
            else
                return std::atoi(env);
        }
    }

    template <typename T, typename... Args>
    class Cache
    {
    private:
        std::function<T(Args...)> func;
        std::size_t size;
        std::deque<std::pair<std::tuple<typename std::remove_reference<Args>::type...>, T>> data;

    public:
        Cache(T (*f)(Args...), std::size_t size) : func(f), size(size){};

        T operator()(Args... args)
        {
            if (this->size == 0)
                return func(args...);

            auto temp = std::tie(args...);

            auto i = std::find_if(
                this->data.cbegin(), this->data.cend(),
                [&](auto &i)
                {
                    return i.first == temp;
                });

            if (i == this->data.cend())
            {
                auto res = std::make_pair(temp, func(args...));

                if (this->data.size() >= this->size)
                    this->data.pop_front();

                return this->data.emplace_back(res).second;
            }

            return i->second;
        }
    };
}
