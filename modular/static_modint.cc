#include <cstdint>
#include <iostream>
#include <type_traits>

namespace internal {

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

struct ModIntBase {};

template <typename T>
using is_modint = std::is_base_of<ModIntBase, T>;

template <typename T>
using is_modint_t = typename std::is_base_of<ModIntBase, T>::type;

}  // namespace internal

template <int kMod>
class LazyMontgomeryModInt : public internal::ModIntBase {
  static_assert(kMod > 0, "Mod must be positive");
  static_assert(kMod % 2 == 1, "Mod must be odd");
  static_assert(kMod < (1 << 30), "Mod must be less than 2^30");

 public:
  using Mint = LazyMontgomeryModInt<kMod>;

  constexpr LazyMontgomeryModInt() : value_(0) {}

  template <typename T,
            typename std::enable_if_t<std::is_signed_v<T>>* = nullptr>
  constexpr LazyMontgomeryModInt(T n)  // NOLINT(*explicit-constructor*)
      : value_(0) {
    if (n %= kMod; n < 0) {
      n += kMod;
    }
    value_ = space_.Transform(static_cast<uint32_t>(n));
  }

  template <typename T,
            typename std::enable_if_t<std::is_unsigned_v<T>>* = nullptr>
  constexpr LazyMontgomeryModInt(T n)  // NOLINT(*explicit-constructor*)
      : value_(space_.Transform(static_cast<uint32_t>(n % kMod))) {}

  static constexpr int Mod() { return kMod; }

  static constexpr Mint Raw(uint32_t value) {
    Mint result;
    result.value_ = space_.Transform(value);
    return result;
  }

  [[nodiscard]] constexpr Mint operator+() const noexcept { return *this; }

  [[nodiscard]] constexpr Mint operator-() const noexcept {
    return Mint() - *this;
  }

  template <typename T,
            typename std::enable_if_t<std::is_integral_v<T>>* = nullptr>
  [[nodiscard]] constexpr explicit operator T() const {
    uint32_t q = space_.Reduce(value_);
    return static_cast<T>(q < kUMod_ ? q : q - kUMod_);
  }

  [[nodiscard]] constexpr Mint Inverse() const noexcept {
    return Power(*this, kUMod_ - 2);
  }

  constexpr Mint& operator+=(const Mint& other) noexcept {
    if (static_cast<int>(value_ += other.value_ - kDoubleUMod_) < 0) {
      value_ += kDoubleUMod_;
    }
    return *this;
  }

  constexpr Mint& operator-=(const Mint& other) noexcept {
    if (static_cast<int>(value_ -= other.value_) < 0) {
      value_ += kDoubleUMod_;
    }
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
    return Mint(lhs) += rhs;
  }

  constexpr friend Mint operator-(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) -= rhs;
  }

  constexpr friend Mint operator*(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) *= rhs;
  }

  constexpr friend Mint operator/(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) /= rhs;
  }

  constexpr friend bool operator==(const Mint& lhs, const Mint& rhs) {
    return (lhs.value_ < kUMod_ ? lhs.value_ : lhs.value_ - kUMod_) ==
           (rhs.value_ < kUMod_ ? rhs.value_ : rhs.value_ - kUMod_);
  }

  constexpr friend bool operator!=(const Mint& lhs, const Mint& rhs) {
    return (lhs.value_ < kUMod_ ? lhs.value_ : lhs.value_ - kUMod_) !=
           (rhs.value_ < kUMod_ ? rhs.value_ : rhs.value_ - kUMod_);
  }

  constexpr friend Mint Power(Mint mint, uint64_t n) noexcept {
    Mint result = Raw(1U);
    while (n > 0) {
      if (n % 2 == 1) {
        result *= mint;
      }
      mint *= mint;
      n /= 2;
    }
    return result;
  }

  friend std::istream& operator>>(std::istream& istream, Mint& mint) {
    uint32_t value;
    istream >> value;
    mint = Mint::Raw(value);
    return istream;
  }

  friend std::ostream& operator<<(std::ostream& ostream, const Mint& mint) {
    return ostream << static_cast<int>(mint);
  }

 private:
  static constexpr uint32_t kUMod_ = kMod;
  static constexpr uint32_t kDoubleUMod_ = kUMod_ << 1;
  static constexpr internal::Montgomery space_ = internal::Montgomery(kUMod_);

  uint32_t value_;
};
