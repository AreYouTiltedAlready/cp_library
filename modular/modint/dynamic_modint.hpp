// #include "modular/modint/modint_base.hpp"

namespace modular {
namespace modint {
namespace dynamic_modint {

template <typename T, int kId>
class DynamicModInt : public internal::ModIntBase {
 public:
  using signed_t = std::make_signed_t<T>;
  using unsigned_t = std::make_unsigned_t<T>;
  static void SetMod(T new_mod) {
    unsigned_mod_ = new_mod;
    barrett_ = barrett::Barrett<unsigned_t>(new_mod);
  }

  using Mint = DynamicModInt<T, kId>;

  constexpr DynamicModInt() : value_(0) {}

  template <typename U>
  requires(!unsigned_int_or_int64<U>)
      DynamicModInt(U n)  // NOLINT(*explicit-constructor*)
      : value_(0) {
    if (n %= unsigned_mod_; n < 0) {
      n += unsigned_mod_;
    }
    value_ = static_cast<unsigned_t>(n);
  }

  template <unsigned_int_or_int64 U>
  DynamicModInt(U n)  // NOLINT(*explicit-constructor*)
      : value_(static_cast<unsigned_t>(n % unsigned_mod_)) {}

  static unsigned_t UMod() { return unsigned_mod_; }

  static Mint Raw(unsigned_t value) {
    Mint result;
    result.value_ = value;
    return result;
  }

  [[nodiscard]] Mint operator+() const noexcept { return *this; }

  [[nodiscard]] Mint operator-() const noexcept { return Mint() - *this; }

  [[nodiscard]] unsigned_t Get() const { return value_; }

  template <std::integral U>
  [[nodiscard]] explicit operator U() const {
    return static_cast<U>(value_);
  }

  [[nodiscard]] Mint Inverse() const noexcept {
    return Power(*this, unsigned_mod_ - 2);
  }

  Mint& operator+=(const Mint& other) noexcept {
    if ((value_ += other.value_) >= unsigned_mod_) {
      value_ -= unsigned_mod_;
    }
    return *this;
  }

  Mint& operator-=(const Mint& other) noexcept {
    if ((value_ += unsigned_mod_ - other.value_) >= unsigned_mod_) {
      value_ -= unsigned_mod_;
    }
    return *this;
  }

  Mint& operator*=(const Mint& other) noexcept {
    value_ = barrett_.Product(value_, other.value_);
    return *this;
  }

  Mint& operator/=(const Mint& other) noexcept {
    value_ = barrett_.Product(value_, other.Inverse().value_);
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
  static inline unsigned_t unsigned_mod_ = 1U;
  static inline barrett::Barrett<unsigned_t> barrett_ =
      barrett::Barrett<unsigned_t>(unsigned_mod_);

  unsigned_t value_;
};

}  // namespace dynamic_modint

template <int kMod>
using DynamicMInt = dynamic_modint::StaticModInt<int, kMod>;

template <int64_t kMod>
using DynamicMLong = dynamic_modint::StaticModInt<int64_t, kMod>;

}  // namespace modint
}  // namespace modular
