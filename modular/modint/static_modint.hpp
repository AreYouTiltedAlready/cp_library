// #include "modular/modint/modint_base.hpp"

namespace modular {
namespace modint {
namespace static_modint {

template <typename T, T kMod>
class StaticModInt : public internal::ModIntBase {
 public:
  using signed_t = std::make_signed_t<T>;
  using unsigned_t = std::make_unsigned_t<T>;
  using Mint = StaticModInt<T, kMod>;

  constexpr StaticModInt() : value_(0) {}

  template <typename U>
  requires(!unsigned_int_or_int64<U>) constexpr StaticModInt(
      U n)  // NOLINT(*explicit-constructor*)
      : value_(0) {
    if (n %= kMod; n < 0) {
      n += kMod;
    }
    value_ = static_cast<unsigned_t>(n);
  }

  template <unsigned_int_or_int64 U>
  constexpr StaticModInt(U n)  // NOLINT(*explicit-constructor*)
      : value_(static_cast<unsigned_t>(n % kMod)) {}

  static constexpr unsigned_t UMod() { return kUnsignedMod; }

  static constexpr Mint Raw(unsigned_t value) {
    Mint result;
    result.value_ = value;
    return result;
  }

  [[nodiscard]] constexpr Mint operator+() const noexcept { return *this; }

  [[nodiscard]] constexpr Mint operator-() const noexcept {
    return Mint() - *this;
  }

  [[nodiscard]] constexpr unsigned_t Get() const { return value_; }

  template <std::integral U>
  [[nodiscard]] constexpr explicit operator U() const {
    return static_cast<U>(value_);
  }

  [[nodiscard]] constexpr Mint Inverse() const noexcept {
    return Power(*this, kUnsignedMod - 2);
  }

  constexpr Mint& operator+=(const Mint& other) noexcept {
    if ((value_ += other.value_) >= kUnsignedMod) {
      value_ -= kUnsignedMod;
    }
    return *this;
  }

  constexpr Mint& operator-=(const Mint& other) noexcept {
    if ((value_ += kMod - other.value_) >= kUnsignedMod) {
      value_ -= kUnsignedMod;
    }
    return *this;
  }

  constexpr Mint& operator*=(const Mint& other) noexcept {
    value_ = kBarrett.Product(value_, other.value_);
    return *this;
  }

  constexpr Mint& operator/=(const Mint& other) noexcept {
    value_ = kBarrett.Product(value_, other.Inverse().value_);
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
    return lhs.value_ == rhs.value_;
  }

  constexpr friend bool operator!=(const Mint& lhs, const Mint& rhs) {
    return lhs.value_ != rhs.value_;
  }

  constexpr friend Mint Power(Mint mint, uint64_t n) noexcept {
    Mint res = Raw(1U);
    while (n > 0) {
      if (n % 2 == 1) {
        res *= mint;
      }
      mint *= mint;
      n /= 2;
    }
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
  static constexpr unsigned_t kUnsignedMod = kMod;
  static constexpr barrett::Barrett<unsigned_t> kBarrett =
      barrett::Barrett<unsigned_t>(kUnsignedMod);

  unsigned_t value_;
};

}  // namespace static_modint

template <int kMod>
using StaticMInt = static_modint::StaticModInt<int, kMod>;

template <int64_t kMod>
using StaticMLong = static_modint::StaticModInt<int64_t, kMod>;


}  // namespace modint
}  // namespace modular
