#include <cstdint>

class Montgomery {
 public:
  constexpr explicit Montgomery(uint32_t mod)
      : mod_(mod),
        r_square_((static_cast<uint64_t>(1) << 32) % mod),
        mod_inverse_(1) {
    r_square_ = static_cast<uint64_t>(r_square_) * r_square_ % mod;
    for (int i = 0; i < 5; ++i) {
      mod_inverse_ *= static_cast<uint32_t>(2) - mod * mod_inverse_;
    }
  }

  [[nodiscard]] constexpr uint32_t Product(uint32_t lhs, uint32_t rhs) const {
    return Reduce(static_cast<uint64_t>(lhs) * rhs);
  }

  // returns a value congruent to mod in range [0, 2 * mod)
  [[nodiscard]] constexpr uint32_t Reduce(uint64_t n) const {
    uint32_t q = static_cast<uint32_t>(n) * mod_inverse_;
    uint32_t m = (static_cast<uint64_t>(q) * mod_) >> 32;
    return (n >> 32) + mod_ - m;
  }

  [[nodiscard]] constexpr uint32_t Transform(uint32_t n) const {
    return Product(n, r_square_);
  }

 private:
  uint32_t mod_;
  uint32_t r_square_;
  uint32_t mod_inverse_;
};
