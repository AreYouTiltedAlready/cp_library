#include <bits/stdc++.h>

namespace {

template <class Fun>
class y_combinator_result {
  Fun fun_;

 public:
  template <class T>
  explicit y_combinator_result(T&& fun) : fun_(std::forward<T>(fun)) {}

  template <class... Args>
  decltype(auto) operator()(Args&&... args) {
    return fun_(std::ref(*this), std::forward<Args>(args)...);
  }
};

template <class Fun>
decltype(auto) y_combinator(Fun&& fun) {
  return y_combinator_result<std::decay_t<Fun>>(std::forward<Fun>(fun));
}

#include <type_traits>
#include <vector>

template <typename T>
void FWHT(std::vector<T>& v) {
  const int n = static_cast<int>(v.size());
  for (int i = 1; i < n; i *= 2) {
    for (int j = 0; j < n; j += i * 2) {
      for (int k = 0; k < i; ++k) {
        T x = v[j + k];
        T y = v[i + j + k];
        v[j + k] = x + y;
        v[i + j + k] = x - y;
      }
    }
  }
}

template <typename T, typename U>
std::vector<std::common_type_t<T, U>> XorConvolution(std::vector<T> lhs,
                                                     std::vector<U> rhs) {
  const int n = static_cast<int>(lhs.size());
  FWHT(lhs);
  FWHT(rhs);
  using V = std::common_type_t<T, U>;
  std::vector<V> result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = static_cast<V>(lhs[i]) * rhs[i];
  }
  FWHT(result);
  for (V& c : result) {
    c /= n;
  }
  return result;
}

class Barrett {
 public:
  constexpr explicit Barrett(uint32_t mod)
      : mod_inverse(static_cast<uint64_t>(-1) / mod + 1), mod(mod) {}

  [[nodiscard]] uint32_t Product(uint32_t lhs, uint32_t rhs) const {
    return (*this)(static_cast<uint64_t>(lhs) * rhs);
  }

  [[nodiscard]] uint32_t operator()(uint64_t n) const {
    auto x = static_cast<uint64_t>(
        (static_cast<__uint128_t>(n) * mod_inverse) >> 64);
    uint64_t m = x * mod;
    return n - m + (n < m ? mod : 0);
  }

 private:
  uint64_t mod_inverse;
  uint32_t mod;
};

// Static modint class (for 32-bit compile-time modulos)
// Uses barrett reduction for multiplication. Quite fast in practice.
// Usage:
// using mint = StaticModint<998244353>; // whatever (1000000007, ...)
template <uint32_t kMod>
class StaticModint {
 public:
  static_assert(kMod < (1U << 30));

  StaticModint() : value_(0) {}

  template <typename T,
            typename std::enable_if_t<
                std::is_integral_v<T> && std::is_signed_v<T>, void>* = nullptr>
  StaticModint(T value)
      : value_(  // NOLINT(*-explicit-constructor)
            value % kMod) {
    if (value < 0) {
      value += kMod;
    }
  }

  template <typename T, typename std::enable_if_t<std::is_integral_v<T> &&
                                                      std::is_unsigned_v<T>,
                                                  void>* = nullptr>
  StaticModint(T value)
      : value_(  // NOLINT(*-explicit-constructor)
            value % kMod) {}

  using Mint = StaticModint<kMod>;

  static Mint Raw(unsigned value) {
    Mint result;
    result.value_ = value;
    return result;
  }

  [[nodiscard]] Mint operator+() const noexcept { return Mint(*this); }

  [[nodiscard]] Mint operator-() const noexcept { return Mint(kMod - value_); }

  [[nodiscard]] bool IsZero() const noexcept { return value_ == 0; }

  [[nodiscard]] explicit operator uint32_t() const { return value_; }

  [[nodiscard]] Mint Inverse() const noexcept {
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
      if (n % 2 == 1) {
        result *= m;
      }
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

// https://judge.yosupo.jp/problem/bitwise_xor_convolution
// https://judge.yosupo.jp/submission/179442
void RunCase([[maybe_unused]] int testcase) {
  int n;
  std::cin >> n;
  const int m = 1 << n;

  using Mint = StaticModint<998244353>;
  std::vector<Mint> p(m);
  for (Mint& z : p) {
    std::cin >> z;
  }
  std::vector<Mint> q(m);
  for (Mint& z : q) {
    std::cin >> z;
  }

  std::vector<Mint> r = XorConvolution(p, q);
  for (const Mint& z : r) {
    std::cout << z << " ";
  }
}

void Main() {
  int testcases = 1;
  // std::cin >> testcases;
  for (int tt = 1; tt <= testcases; ++tt) {
    RunCase(tt);
  }
}

}  // namespace

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);
  Main();
  return 0;
}
