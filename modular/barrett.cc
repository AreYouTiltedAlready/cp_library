#include <cstdint>
#include <limits>

class Barrett {
public:
  constexpr explicit Barrett(uint32_t mod)
      : mod_inverse(static_cast<uint64_t>(-1) / mod + 1), mod(mod) {}

  [[nodiscard]]
  uint32_t Product(uint32_t lhs, uint32_t rhs) const {
    return (*this)(static_cast<uint64_t>(lhs) * rhs);
  }

  [[nodiscard]]
  uint32_t operator()(uint64_t n) const {
    auto x = static_cast<uint64_t>(
        (static_cast<__uint128_t>(n) * mod_inverse)
            >> 64);
    uint64_t m = x * mod;
    return n - m + (n < m ? mod : 0);
  }

private:
  uint64_t mod_inverse;
  uint32_t mod;
};

using uint128_t = __uint128_t;

class Barrett128 {
public:
  constexpr explicit Barrett128(uint64_t mod)
      : mod(mod), mod_inverse_high((static_cast<uint128_t>(-1) / mod + 1) >> 64)
        , mod_inverse_low((static_cast<uint128_t>(-1) / mod + 1) -
                          (static_cast<uint128_t>(mod_inverse_high) << 64)) {
  }

  [[nodiscard]]
  uint64_t Product(uint64_t lhs, uint64_t rhs) const {
    return (*this)(static_cast<uint128_t>(lhs) * rhs);
  }

  [[nodiscard]]
  uint64_t operator()(uint128_t n) const {
    uint128_t x = Product(n & std::numeric_limits<uint64_t>::max(), n >> 64,
                          mod_inverse_low,
                          mod_inverse_high);
    uint128_t m = x * mod;
    return n - m + (n < m ? mod : 0);
  }

private:
  static uint128_t
  Product(uint64_t lhs_low, uint64_t lhs_high, uint64_t rhs_low,
          uint64_t rhs_high) {
    // (a*2^64 + b) * (c*2^64 + d) =
    // (a*c) * 2^128 + (a*d + b*c)*2^64 + (b*d)
    uint128_t ac = static_cast<uint128_t>(lhs_high) * rhs_high;
    uint128_t ad = static_cast<uint128_t>(lhs_high) * rhs_low;
    uint128_t bc = static_cast<uint128_t>(lhs_low) * rhs_high;
    uint128_t bd = static_cast<uint128_t>(lhs_low) * rhs_low;
    uint128_t carry = static_cast<uint128_t>(static_cast<uint64_t>(ad)) +
                      static_cast<uint128_t>(static_cast<uint64_t>(bc) +
                                             (bd >> 64));
    return ac + (ad >> 64) + (bc >> 64) + (carry >> 64);
  }

  uint64_t mod;
  uint64_t mod_inverse_high;
  uint64_t mod_inverse_low;
};
