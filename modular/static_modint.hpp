#include <cstdint>
#include <iostream>
#include <limits>
#include <type_traits>

namespace modular {

namespace internal {
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
        mod_inverse_(1) {
    r_square_ = static_cast<promoted_t>(r_square_) * r_square_ % mod;
    for (int i = 0; i < 6; ++i) {
      mod_inverse_ *= static_cast<unsigned_t>(2) - mod_ * mod_inverse_;
    }
  }

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

  [[nodiscard]] constexpr bool AreEqual(unsigned_t lhs, unsigned_t rhs) const {
    return (lhs < mod_ ? lhs : lhs - mod_) == (rhs < mod_ ? rhs : rhs - mod_);
  }

  [[nodiscard]] constexpr unsigned_t Power(unsigned_t x, uint64_t n) const {
    unsigned_t result = Transform(1U);
    while (n > 0) {
      if (n % 2 == 1) {
        result = Product(result, x);
      }
      x = Product(x, x);
      n /= 2;
    }
    return result;
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
};

}  // namespace montgomery

struct ModIntBase {};

template <typename T>
using is_modint = std::is_base_of<ModIntBase, T>;

template <typename T>
using is_modint_t = typename std::is_base_of<ModIntBase, T>::type;

template <typename T, T kMod>
class LazyMontgomeryModInt : public internal::ModIntBase {
  using signed_t = std::make_signed_t<T>;
  using unsigned_t = std::make_unsigned_t<T>;

  static_assert(kMod > 0, "Mod must be positive");
  static_assert(kMod % 2 == 1, "Mod must be odd");
  static_assert((kMod < (static_cast<uint64_t>(1)
                         << (std::numeric_limits<signed_t>::digits - 1))),
                "Mod is too large");

 public:
  using Mint = LazyMontgomeryModInt<T, kMod>;

  constexpr LazyMontgomeryModInt() : value_(0) {}

  template <typename U,
            typename std::enable_if_t<std::is_signed_v<U>>* = nullptr>
  constexpr LazyMontgomeryModInt(U n)  // NOLINT(*explicit-constructor*)
      : value_(0) {
    if (n %= kMod; n < 0) {
      n += kMod;
    }
    value_ = space_.Transform(static_cast<unsigned_t>(n));
  }

  template <typename U,
            typename std::enable_if_t<std::is_unsigned_v<U>>* = nullptr>
  constexpr LazyMontgomeryModInt(U n)  // NOLINT(*explicit-constructor*)
      : value_(space_.Transform(static_cast<unsigned_t>(n % kMod))) {}

  static constexpr unsigned_t UMod() { return kUMod_; }

  static constexpr Mint Raw(unsigned_t value) {
    Mint result;
    result.value_ = space_.Transform(value);
    return result;
  }

  [[nodiscard]] constexpr Mint operator+() const noexcept { return *this; }

  [[nodiscard]] constexpr Mint operator-() const noexcept {
    return Mint() - *this;
  }

  template <typename U,
            typename std::enable_if_t<std::is_integral_v<U>>* = nullptr>
  [[nodiscard]] constexpr explicit operator U() const {
    unsigned_t q = space_.Reduce(value_);
    return static_cast<T>(q < kUMod_ ? q : q - kUMod_);
  }

  [[nodiscard]] constexpr Mint Inverse() const noexcept {
    return Power(*this, kUMod_ - 2);
  }

  constexpr Mint& operator+=(const Mint& other) noexcept {
    value_ = space_.Sum(value_, other.value_);
    return *this;
  }

  constexpr Mint& operator-=(const Mint& other) noexcept {
    value_ = space_.Difference(value_, other.value_);
    return *this;
  }

  constexpr Mint& operator*=(const Mint& other) noexcept {
    value_ = space_.Product(value_, other.value_);
    return *this;
  }

  constexpr Mint& operator/=(const Mint& other) noexcept {
    value_ = space_.Product(value_, other.Inverse().value_);
    return *this;
  }

  constexpr friend Mint operator+(const Mint& lhs, const Mint& rhs) {
    Mint res{};
    res.value_ = space_.Sum(lhs.value_, rhs.value_);
    return res;
  }

  constexpr friend Mint operator-(const Mint& lhs, const Mint& rhs) {
    Mint res{};
    res.value_ = space_.Difference(lhs.value_, rhs.value_);
    return res;
  }

  constexpr friend Mint operator*(const Mint& lhs, const Mint& rhs) {
    Mint res{};
    res.value_ = space_.Product(lhs.value_, rhs.value_);
    return res;
  }

  constexpr friend Mint operator/(const Mint& lhs, const Mint& rhs) {
    Mint res{};
    res.value_ = space_.Product(lhs.value_, rhs.Inverse().value_);
    return res;
  }

  constexpr friend bool operator==(const Mint& lhs, const Mint& rhs) {
    return space_.AreEqual(lhs.value_, rhs.value_);
  }

  constexpr friend bool operator!=(const Mint& lhs, const Mint& rhs) {
    return !space_.AreEqual(lhs.value_, rhs.value_);
  }

  constexpr friend Mint Power(Mint mint, uint64_t n) noexcept {
    Mint res;
    res.value_ = space_.Power(mint.value_, n);
    return res;
  }

  friend std::istream& operator>>(std::istream& istream, Mint& mint) {
    unsigned_t value;
    istream >> value;
    mint = Mint::Raw(value);
    return istream;
  }

  friend std::ostream& operator<<(std::ostream& ostream, const Mint& mint) {
    return ostream << static_cast<signed_t>(mint);
  }

 private:
  static constexpr unsigned_t kUMod_ = kMod;
  static constexpr internal::montgomery::MontgomerySpace<unsigned_t> space_ =
      internal::montgomery::MontgomerySpace<unsigned_t>(kUMod_);

  unsigned_t value_;
};

}  // namespace internal

template <int kMod>
using MInt = internal::LazyMontgomeryModInt<int, kMod>;

template <int64_t kMod>
using MLong = internal::LazyMontgomeryModInt<int64_t, kMod>;

}  // namespace modular
