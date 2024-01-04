#include <cstdint>
#include <type_traits>

namespace modular {
namespace integral_type_traits {

template <typename T>
struct integral_promotion {
  using type = T;
};

template <>
struct integral_promotion<int> {
  using type = int64_t;
};

template <>
struct integral_promotion<uint32_t> {
  using type = uint64_t;
};

template <>
struct integral_promotion<int64_t> {
  using type = __int128_t;
};

template <>
struct integral_promotion<uint64_t> {
  using type = __uint128_t;
};

}  // namespace integral_type_traits

template <typename T>
using integral_promotion_t =
    typename integral_type_traits::integral_promotion<T>::type;

template <typename T>
concept unsigned_int_or_int64 =
    std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>;

}  // namespace modular
