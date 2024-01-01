#include <bits/stdc++.h>

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

namespace bit {

uint32_t bit_ceil(uint32_t n) {
  uint32_t res = 1;
  while (res < n) {
    res *= 2;
  }
  return res;
}

uint32_t countr_zero(uint32_t n) { return __builtin_ctz(n); }
uint32_t countl_zero(uint32_t n) { return __builtin_clz(n); }

constexpr uint32_t countl_zero_constexpr(uint32_t n) {
  if (n == 0) {
    return 32U;
  }
  uint32_t res = 31U;
  while ((n >> res) % 2 == 0) {
    res -= 1;
  }
  return 31U - res;
}

constexpr uint32_t countr_zero_constexpr(uint32_t n) {
  uint32_t res = 0;
  while ((n >> res) % 2 == 0) {
    res += 1;
  }
  return res;
}

}  // namespace bit

constexpr int PowModConstexpr(int x, int n, int mod) {
  int result = 1;
  while (n > 0) {
    if (n % 2 == 1) {
      result = static_cast<int>(static_cast<int64_t>(result) * x % mod);
    }
    x = static_cast<int>(static_cast<int64_t>(x) * x % mod);
    n /= 2;
  }
  return result;
}

// mod must be prime
constexpr int PrimitiveRootConstexpr(int mod) {
  if (mod == 2) {
    return 1;
  }
  if (mod == 167772161) {
    return 3;
  }
  if (mod == 469762049) {
    return 3;
  }
  if (mod == 754974721) {
    return 11;
  }
  if (mod == 998244353) {
    return 3;
  }

  int divisors[20] = {};
  divisors[0] = 2;
  int count = 1;
  int d =
      (mod - 1) >> bit::countr_zero_constexpr(static_cast<uint32_t>(mod - 1));

  for (int i = 3; static_cast<int64_t>(i) * i <= d; i += 2) {
    if (d % i == 0) {
      divisors[count++] = i;
      while (d % i == 0) {
        d /= i;
      }
    }
  }
  if (d > 1) {
    divisors[count++] = d;
  }

  for (int g = 2;; g++) {
    bool ok = true;
    for (int i = 0; i < count; i++) {
      if (PowModConstexpr(g, (mod - 1) / divisors[i], mod) == 1) {
        ok = false;
        break;
      }
    }
    if (ok) {
      return g;
    }
  }
}

template <int kMod>
constexpr int kPrimitiveRoot = PrimitiveRootConstexpr(kMod);

template <int kMaxRank>
class NTT {
 public:
  explicit NTT() : bit_inverse_() {
    for (int i = 1; i < 1 << kMaxRank; ++i) {
      bit_inverse_[i] =
          (bit_inverse_[i >> 1] >> 1) + ((i % 2) << (kMaxRank - 1));
    }
  }

  template <int kMod = 998244353, class T,
            std::enable_if_t<std::is_integral<T>::value>* = nullptr>
  std::vector<T> Convolution(const std::vector<T>& lhs,
                             const std::vector<T>& rhs) {
    const int n = static_cast<int>(lhs.size());
    const int m = static_cast<int>(rhs.size());
    if (std::min(n, m) == 0) {
      return {};
    }

    using Mint = LazyMontgomeryModInt<kMod>;

    const int z =
        static_cast<int>(bit::bit_ceil(static_cast<uint32_t>((n + m - 1))));
    assert((Mint::Mod() - 1) % z == 0);

    std::vector<Mint> lhs_copy(n);
    std::vector<Mint> rhs_copy(m);
    std::copy(lhs.cbegin(), lhs.cend(), lhs_copy.begin());
    std::copy(rhs.cbegin(), rhs.cend(), rhs_copy.begin());

    auto c = Convolution(std::move(lhs_copy), std::move(rhs_copy));
    std::vector<T> result(n + m - 1);
    for (int i = 0; i < n + m - 1; i++) {
      result[i] = static_cast<T>(c[i]);
    }
    return result;
  }

  template <class Mint, internal::is_modint_t<Mint>* = nullptr>
  std::vector<Mint> Convolution(const std::vector<Mint>& lhs,
                                const std::vector<Mint>& rhs) {
    const int n = static_cast<int>(lhs.size());
    const int m = static_cast<int>(rhs.size());
    if (std::min(n, m) == 0) {
      return {};
    }

    const int z =
        static_cast<int>(bit::bit_ceil(static_cast<uint32_t>((n + m - 1))));
    assert((Mint::Mod() - 1) % z == 0);

    if (std::min(n, m) <= 60) {
      return ConvolutionNaive(lhs, rhs);
    }

    return ConvolutionNTT(lhs, rhs);
  }

  template <class Mint, internal::is_modint_t<Mint>* = nullptr>
  std::vector<Mint> Convolution(std::vector<Mint>&& lhs,
                                std::vector<Mint>&& rhs) {
    const int n = static_cast<int>(lhs.size());
    const int m = static_cast<int>(rhs.size());
    if (std::min(n, m) == 0) {
      return {};
    }

    const int z =
        static_cast<int>(bit::bit_ceil(static_cast<uint32_t>((n + m - 1))));
    assert((Mint::Mod() - 1) % z == 0);

    if (std::min(n, m) <= 60) {
      return ConvolutionNaive(lhs, rhs);
    }

    return ConvolutionNTT(lhs, rhs);
  }

 private:
  template <class Mint, int kG = kPrimitiveRoot<Mint::Mod()>,
            internal::is_modint_t<Mint>* = nullptr>
  struct NTTInfo {
    static constexpr int kBinaryRank =
        static_cast<int>(bit::countr_zero_constexpr(Mint::Mod() - 1));

    constexpr NTTInfo() : roots() {
      roots[kBinaryRank - 1] =
          Power(Mint(kG), (Mint::Mod() - 1) >> kBinaryRank);
      for (int i = kBinaryRank - 2; i >= 0; --i) {
        roots[i] = roots[i + 1] * roots[i + 1];
      }
    }

    std::array<Mint, kBinaryRank> roots;
  };

  template <class Mint, internal::is_modint_t<Mint>* = nullptr>
  std::vector<Mint> ConvolutionNTT(std::vector<Mint> lhs,
                                   std::vector<Mint> rhs) {
    const int n = static_cast<int>(lhs.size());
    const int m = static_cast<int>(rhs.size());
    const int z =
        static_cast<int>(bit::bit_ceil(static_cast<uint32_t>((n + m - 1))));
    lhs.resize(z);
    Butterfly(lhs);
    rhs.resize(z);
    Butterfly(rhs);
    for (int i = 0; i < z; i++) {
      lhs[i] *= rhs[i];
    }
    Butterfly(lhs);
    std::reverse(lhs.begin() + 1, lhs.end());
    lhs.resize(n + m - 1);
    Mint z_inverse = Mint(z).Inverse();
    for (Mint& c : lhs) {
      c *= z_inverse;
    }
    return lhs;
  }

  template <class Mint, internal::is_modint_t<Mint>* = nullptr>
  std::vector<Mint> ConvolutionNaive(const std::vector<Mint>& lhs,
                                     const std::vector<Mint>& rhs) {
    const int n = static_cast<int>(lhs.size());
    const int m = static_cast<int>(rhs.size());
    std::vector<Mint> res(n + m - 1);
    if (n < m) {
      for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
          res[i + j] += lhs[i] * rhs[j];
        }
      }
    } else {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
          res[i + j] += lhs[i] * rhs[j];
        }
      }
    }
    return res;
  }

  template <typename Mint, internal::is_modint_t<Mint>* = nullptr>
  void Butterfly(std::vector<Mint>& data) {
    const int n = static_cast<int>(data.size());
    const int rank =
        static_cast<int>(bit::countr_zero(static_cast<uint32_t>(n)));

    {
      const int shift = kMaxRank - rank;
      for (int i = 0; i < n; ++i) {
        if (int inv_i = bit_inverse_[i] >> shift; inv_i < i) {
          std::swap(data[i], data[inv_i]);
        }
      }
    }

    static constexpr NTTInfo<Mint> kNTTInfo{};
    for (int i = 0; i < rank; ++i) {
      const int length = 1 << i;
      const Mint root = kNTTInfo.roots[i];
      for (int base = 0; base < n; base += length * 2) {
        Mint w = Mint::Raw(1U);
        for (int offset = 0; offset < length; ++offset) {
          Mint z = data[base + offset + length] * w;
          data[base + offset + length] = data[base + offset] - z;
          data[base + offset] = data[base + offset] + z;
          w *= root;
        }
      }
    }
  }

  std::array<int, kMaxRank> bit_inverse_;
};
