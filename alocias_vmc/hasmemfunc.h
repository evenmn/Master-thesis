#ifndef HASMEMFUNC_H
#define HASMEMFUNC_H

#include <type_traits>
        
// macro for template struct that checks for member function on compile time
// and keeps a boolean 'value' record (true if exists false if not). NOTE: Does
// not work with template functions (that was a different hell to dive into...)
#define HAS_MEM_FUNC(func, name) \
    template <typename T, typename... Args> struct name { \
        private: \
            template<typename C, typename = \
                decltype(std::declval<C>().func(std::declval<Args>()...))> \
            static std::true_type test(int); \
            template<typename C> static std::false_type test(...); \
        public: \
            static constexpr bool value = decltype(test<T>(0))::value; \
    };

// macro which expands above macro for a function func which is checked for and
// set to name and returns the return from func if it exists and returns the
// sendt in object t of type T.
#define MEM_FUNC(func, check, name) \
    HAS_MEM_FUNC(func, check); \
    template<typename T, typename... Args> auto name(T* t, Args... args) { \
        if constexpr (check<T,Args...>::value) { \
            return t->func(args...); \
        } else { \
            return t; \
        } \
    }

#endif /* HASMEMFUNC_H */
