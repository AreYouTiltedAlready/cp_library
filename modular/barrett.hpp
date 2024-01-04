// #include "modular/type_traits.hpp"

namespace modular {
namespace barrett {

using uint128_t = __uint128_t;

constexpr uint128_t MultiplyHigh(uint128_t lhs, uint128_t rhs) {
  const auto a = static_cast<uint64_t>(lhs >> 64);
  const auto b = static_cast<uint64_t>(lhs);
  const auto c = static_cast<uint64_t>(rhs >> 64);
  const auto d = static_cast<uint64_t>(rhs);
  uint128_t ac = static_cast<uint128_t>(a) * c;
  uint128_t ad = static_cast<uint128_t>(a) * d;
  uint128_t bc = static_cast<uint128_t>(b) * c;
  uint128_t bd = static_cast<uint128_t>(b) * d;
  uint128_t carry = static_cast<uint128_t>(static_cast<uint64_t>(ad)) +
                    static_cast<uint128_t>(static_cast<uint64_t>(bc)) +
                    (bd >> 64);
  return ac + (ad >> 64) + (bc >> 64) + (carry >> 64);
}

template <unsigned_int_or_int64 T>
class Barrett {
 public:
  using promoted_t = integral_promotion_t<T>;

  constexpr explicit Barrett(T mod)
      : mod_inverse_(static_cast<promoted_t>(-1) / mod + 1), mod_(mod) {}

  [[nodiscard]] constexpr T Product(T lhs, T rhs) const {
    return Reduce(static_cast<promoted_t>(lhs) * rhs);
  }

  [[nodiscard]] constexpr T Reduce(promoted_t n) const {
    promoted_t m = mod_;
    if constexpr (std::is_same_v<T, uint32_t>) {
      m *= static_cast<promoted_t>((static_cast<uint128_t>(n) * mod_inverse_) >>
                                   64);
    } else {
      m *= MultiplyHigh(n, mod_inverse_);
    }
    n += (n < m ? mod_ : 0);
    return n - m;
  }

 private:
  promoted_t mod_inverse_;
  T mod_;
};

}  // namespace barrett
}  // namespace modular
