// #include "modular/type_traits.hpp"

namespace montgomery {

template <typename T>
class MontgomerySpace {
 public:
  using signed_t = std::make_signed_t<T>;
  using unsigned_t = std::make_unsigned_t<T>;
  using promoted_t = integral_promotion_t<unsigned_t>;
  static constexpr int kTBitWidth = std::numeric_limits<unsigned_t>::digits;

  constexpr explicit MontgomerySpace(T mod)
      : mod_(mod),
        r_square_((static_cast<promoted_t>(1) << kTBitWidth) % mod),
        mod_inverse_(1),
        unit_(),
        neg_unit_() {
    r_square_ = static_cast<promoted_t>(r_square_) * r_square_ % mod;
    for (int i = 0; i < 6; ++i) {
      mod_inverse_ *= static_cast<unsigned_t>(2) - mod_ * mod_inverse_;
    }
    unit_ = Transform(1);
    neg_unit_ = Transform(mod_ - 1);
  }

  [[nodiscard]] constexpr unsigned_t Unit() const { return unit_; }

  [[nodiscard]] constexpr unsigned_t NegUnit() const { return neg_unit_; }

  [[nodiscard]] constexpr unsigned_t Sum(unsigned_t lhs, unsigned_t rhs) const {
    if (static_cast<signed_t>(lhs += rhs - (mod_ << 1)) < 0) {
      lhs += mod_ << 1;
    }
    return lhs;
  }

  [[nodiscard]] constexpr unsigned_t Difference(unsigned_t lhs,
                                                unsigned_t rhs) const {
    if (static_cast<signed_t>(lhs -= rhs) < 0) {
      lhs += mod_ << 1;
    }
    return lhs;
  }

  [[nodiscard]] constexpr unsigned_t Product(unsigned_t lhs,
                                             unsigned_t rhs) const {
    return Reduce(static_cast<promoted_t>(lhs) * rhs);
  }

  [[nodiscard]] constexpr bool TestEquality(unsigned_t lhs,
                                            unsigned_t rhs) const {
    return (lhs < mod_ ? lhs : lhs - mod_) == (rhs < mod_ ? rhs : rhs - mod_);
  }

  [[nodiscard]] constexpr unsigned_t Power(unsigned_t x, uint64_t n) const {
    unsigned_t result = unit_;
    while (n > 0) {
      if (n % 2 == 1) {
        result = Product(result, x);
      }
      x = Product(x, x);
      n /= 2;
    }
    return result;
  }

  [[nodiscard]] constexpr unsigned_t ExactReduce(promoted_t n) const {
    unsigned_t result = Reduce(n);
    return result - (result >= mod_ ? mod_ : 0);
  }

  // returns a value congruent to mod in range [0, 2 * mod)
  [[nodiscard]] constexpr unsigned_t Reduce(promoted_t n) const {
    unsigned_t q = static_cast<unsigned_t>(n) * mod_inverse_;
    unsigned_t m = (static_cast<promoted_t>(q) * mod_) >> kTBitWidth;
    return (n >> kTBitWidth) + mod_ - m;
  }

  [[nodiscard]] constexpr unsigned_t Transform(unsigned_t n) const {
    return Product(n, r_square_);
  }

 private:
  unsigned_t mod_;
  unsigned_t r_square_;
  unsigned_t mod_inverse_;
  unsigned_t unit_;
  unsigned_t neg_unit_;
};

}  // namespace montgomery