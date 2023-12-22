#include <cstdint>
#include <iostream>
#include <type_traits>

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


template <uint32_t kMod>
class StaticModint {
public:
  static_assert(kMod < (1U << 30));

  StaticModint() : value_(0) {}

  template <typename T, typename std::enable_if_t<
      std::is_integral_v<T> && std::is_signed_v<T>, void>* = nullptr>
  StaticModint(T value) : value_( // NOLINT(*-explicit-constructor)
      value % kMod) {
    if (value < 0) { value += kMod; }
  }

  template <typename T, typename std::enable_if_t<
      std::is_integral_v<T> && std::is_unsigned_v<T>, void>* = nullptr>
  StaticModint(T value) : value_( // NOLINT(*-explicit-constructor)
      value % kMod) {}


  using Mint = StaticModint<kMod>;

  static Mint Raw(unsigned value) {
    Mint result;
    result.value_ = value;
    return result;
  }

  [[nodiscard]]
  Mint operator+() const noexcept {
    return Mint(*this);
  }

  [[nodiscard]]
  Mint operator-() const noexcept {
    return Mint(kMod - value_);
  }

  [[nodiscard]]
  bool IsZero() const noexcept {
    return value_ == 0;
  }

  [[nodiscard]]
  explicit operator uint32_t() const {
    return value_;
  }

  [[nodiscard]]
  Mint Inverse() const noexcept {
    Mint m(*this);
    Mint result = Raw(1U);
#pragma GCC unroll(30)
    for (unsigned i = 0; i < 30; ++i) {
      if (((kMod - 2) >> i) % 2 == 1) {
        result *= m;
      }
      m *= m;
    }
    return result;
  }

  Mint& operator+=(const Mint& other) {
    if (value_ += other.value_; value_ >= kMod) {
      value_ -= kMod;
    }
    return *this;
  }

  Mint& operator-=(const Mint& other) {
    if (value_ += kMod - other.value_; value_ >= kMod) {
      value_ -= kMod;
    }
    return *this;
  }

  Mint& operator*=(const Mint& other) {
    value_ = barrett.Product(value_, other.value_);
    return *this;
  }

  Mint& operator/=(const Mint& other) {
    value_ = barrett.Product(value_, other.Inverse().value_);
    return *this;
  }

  friend Mint operator+(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) += rhs;
  }

  friend Mint operator-(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) -= rhs;
  }

  friend Mint operator*(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) *= rhs;
  }

  friend Mint operator/(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) /= rhs;
  }

  friend Mint Power(Mint m, uint64_t n) {
    Mint result = Raw(1U);
    while (n > 0) {
      if (n % 2 == 1) { result *= m; }
      m *= m;
      n /= 2;
    }
    return result;
  }

  friend bool operator==(const Mint& lhs, const Mint& rhs) {
    return lhs.value_ == rhs.value_;
  }

  friend std::istream& operator>>(std::istream& is, Mint& mint) {
    uint32_t value;
    is >> value;
    mint = Mint(value);
    return is;
  }

  friend std::ostream& operator<<(std::ostream& os, const Mint& mint) {
    return os << mint.value_;
  }

private:
  static constexpr Barrett barrett = Barrett(kMod);
  uint32_t value_;
};
