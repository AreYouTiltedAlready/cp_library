#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

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
concept ModInt = std::is_base_of_v<ModIntBase, T>;

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

  friend std::istream& operator>>(std::istream& istream, Mint& Mint) {
    unsigned_t value;
    istream >> value;
    Mint = Mint::Raw(value);
    return istream;
  }

  friend std::ostream& operator<<(std::ostream& ostream, const Mint& Mint) {
    return ostream << static_cast<signed_t>(Mint);
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

namespace ntt {

namespace internal {

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
template <std::unsigned_integral T>
constexpr int PrimitiveRootConstexpr(T mod) {
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
  if (mod == 4179340454199820289) {
    return 3;
  }

  int divisors[20] = {};
  divisors[0] = 2;
  int count = 1;
  int d = static_cast<int>(
      (mod - 1) >>
      std::countr_zero(static_cast<std::make_unsigned_t<T>>(mod - 1)));

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

template <std::unsigned_integral T, T kMod>
constexpr int kPrimitiveRoot = PrimitiveRootConstexpr(kMod);

template <modular::internal::ModInt Mint,
          int kG = kPrimitiveRoot<decltype(Mint::UMod()), Mint::UMod()>>
struct NTTInfo {
  static constexpr int kRank = std::countr_zero(Mint::UMod() - 1);

  constexpr NTTInfo()
      : root(), inv_root(), rate_2(), inv_rate_2(), rate_3(), inv_rate_3() {
    root[kRank] = Power(Mint(kG), (Mint::UMod() - 1) >> kRank);
    inv_root[kRank] = root[kRank].Inverse();
    for (int i = kRank - 1; i >= 0; i--) {
      root[i] = root[i + 1] * root[i + 1];
      inv_root[i] = inv_root[i + 1] * inv_root[i + 1];
    }

    {
      Mint product = 1;
      Mint inv_product = 1;
      for (int i = 0; i <= kRank - 2; i++) {
        rate_2[i] = root[i + 2] * product;
        inv_rate_2[i] = inv_root[i + 2] * inv_product;
        product *= inv_root[i + 2];
        inv_product *= root[i + 2];
      }
    }
    {
      Mint product = 1;
      Mint inv_product = 1;
      for (int i = 0; i <= kRank - 3; i++) {
        rate_3[i] = root[i + 3] * product;
        inv_rate_3[i] = inv_root[i + 3] * inv_product;
        product *= inv_root[i + 3];
        inv_product *= root[i + 3];
      }
    }
  }

  std::array<Mint, kRank + 1> root;      // root[i]^(2^i) == 1
  std::array<Mint, kRank + 1> inv_root;  // root[i] * inv_root[i] == 1

  std::array<Mint, std::max(0, kRank - 2 + 1)> rate_2;
  std::array<Mint, std::max(0, kRank - 2 + 1)> inv_rate_2;

  std::array<Mint, std::max(0, kRank - 3 + 1)> rate_3;
  std::array<Mint, std::max(0, kRank - 3 + 1)> inv_rate_3;
};

template <modular::internal::ModInt Mint>
void Butterfly(std::vector<Mint>& a) {
  const auto n = static_cast<int>(a.size());
  const int rank = std::countr_zero(static_cast<uint32_t>(n));

  static constexpr NTTInfo<Mint> info{};

  int length = 0;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
  while (length < rank) {
    if (rank - length == 1) {
      const int p = 1 << (rank - length - 1);
      Mint root = Mint::Raw(1U);
      for (int s = 0; s < (1 << length); s++) {
        const int offset = s << (rank - length);
        for (int i = 0; i < p; i++) {
          Mint left = a[i + offset];
          Mint right = a[i + offset + p] * root;
          a[i + offset] = left + right;
          a[i + offset + p] = left - right;
        }
        if (s + 1 != (1 << length)) {
          root *= info.rate_2[std::countr_one(static_cast<uint32_t>(s))];
        }
      }
      length += 1;
    } else {
      // 4-base
      const int p = 1 << (rank - length - 2);
      Mint root = Mint::Raw(1U);
      const Mint imag = info.root[2];
      for (int s = 0; s < (1 << length); s++) {
        const Mint root_square = root * root;
        const Mint root_cube = root * root_square;
        const int offset = s << (rank - length);
        for (int i = 0; i < p; i++) {
          const Mint a_0 = a[i + offset];
          const Mint a_1 = a[i + offset + p] * root;
          const Mint a_2 = a[i + offset + 2 * p] * root_square;
          const Mint a_3 = a[i + offset + 3 * p] * root_cube;
          const Mint aux = (a_1 - a_3) * imag;
          a[i + offset] = a_0 + a_2 + a_1 + a_3;
          a[i + offset + 1 * p] = a_0 + a_2 - a_1 - a_3;
          a[i + offset + 2 * p] = a_0 - a_2 + aux;
          a[i + offset + 3 * p] = a_0 - a_2 - aux;
        }
        if (s + 1 != (1 << length)) {
          root *= info.rate_3[std::countr_one(static_cast<uint32_t>(s))];
        }
      }
      length += 2;
    }
  }
}

template <modular::internal::ModInt Mint>
void InverseButterfly(std::vector<Mint>& a) {
  const auto n = static_cast<int>(a.size());
  const int rank = std::countr_zero(static_cast<uint32_t>(n));

  static constexpr NTTInfo<Mint> info{};

  int length = rank;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
  while (length > 0) {
    if (length == 1) {
      const int p = 1 << (rank - length);
      Mint inv_root = Mint::Raw(1U);
      for (int s = 0; s < (1 << (length - 1)); s++) {
        const int offset = s << (rank - length + 1);
        for (int i = 0; i < p; i++) {
          Mint left = a[i + offset];
          Mint right = a[i + offset + p];
          a[i + offset] = left + right;
          a[i + offset + p] = (left - right) * inv_root;
        }
        if (s + 1 != (1 << (length - 1))) {
          inv_root *=
              info.inv_rate_2[std::countr_one(static_cast<uint32_t>(s))];
        }
      }
      length -= 1;
    } else {
      // 4-base
      const int p = 1 << (rank - length);
      Mint inv_root = Mint::Raw(1U);
      const Mint inv_imag = info.inv_root[2];
      for (int s = 0; s < (1 << (length - 2)); s++) {
        const Mint inv_root_square = inv_root * inv_root;
        const Mint inv_root_cube = inv_root * inv_root_square;
        const int offset = s << (rank - length + 2);
        for (int i = 0; i < p; i++) {
          const Mint a_0 = a[i + offset];
          const Mint a_1 = a[i + offset + p];
          const Mint a_2 = a[i + offset + 2 * p];
          const Mint a_3 = a[i + offset + 3 * p];
          const Mint aux = (a_2 - a_3) * inv_imag;
          a[i + offset] = a_0 + a_1 + a_2 + a_3;
          a[i + offset + 1 * p] = (a_0 - a_1 + aux) * inv_root;
          a[i + offset + 2 * p] = (a_0 + a_1 - a_2 - a_3) * inv_root_square;
          a[i + offset + 3 * p] = (a_0 - a_1 - aux) * inv_root_cube;
        }
        if (s + 1 != (1 << (length - 2))) {
          inv_root *=
              info.inv_rate_3[std::countr_one(static_cast<uint32_t>(s))];
        }
      }
      length -= 2;
    }
  }
}

template <modular::internal::ModInt Mint>
std::vector<Mint> ConvolutionNaive(const std::vector<Mint>& lhs,
                                   const std::vector<Mint>& rhs) {
  const auto n = static_cast<int>(lhs.size());
  const auto m = static_cast<int>(rhs.size());
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

template <modular::internal::ModInt Mint>
std::vector<Mint> ConvolutionNTT(std::vector<Mint> lhs, std::vector<Mint> rhs) {
  const auto n = static_cast<int>(lhs.size());
  const auto m = static_cast<int>(rhs.size());
  const auto z =
      static_cast<int>(std::bit_ceil(static_cast<uint32_t>(n + m - 1)));

  lhs.resize(z);
  Butterfly(lhs);
  rhs.resize(z);
  Butterfly(rhs);
  for (int i = 0; i < z; i++) {
    lhs[i] *= rhs[i];
  }
  InverseButterfly(lhs);
  lhs.resize(n + m - 1);
  const Mint z_inverse = Mint(z).Inverse();
  for (Mint& it : lhs) {
    it *= z_inverse;
  }
  return lhs;
}

}  // namespace internal

template <modular::internal::ModInt Mint>
std::vector<Mint> Convolution(std::vector<Mint>&& lhs,
                              std::vector<Mint>&& rhs) {
  const auto n = static_cast<int>(lhs.size());
  const auto m = static_cast<int>(rhs.size());
  if (std::min(n, m) == 0) {
    return {};
  }
  const auto z =
      static_cast<int>(std::bit_ceil(static_cast<uint32_t>(n + m - 1)));
  assert((Mint::UMod() - 1) % z == 0);

  // if (std::min(n, m) <= 60) {
  //   return internal::ConvolutionNaive(lhs, rhs);
  // }
  return internal::ConvolutionNTT(std::move(lhs), std::move(rhs));
}

template <modular::internal::ModInt Mint>
std::vector<Mint> Convolution(const std::vector<Mint>& lhs,
                              const std::vector<Mint>& rhs) {
  const auto n = static_cast<int>(lhs.size());
  const auto m = static_cast<int>(rhs.size());
  if (std::min(n, m) == 0) {
    return {};
  }
  const auto z =
      static_cast<int>(std::bit_ceil(static_cast<uint32_t>(n + m - 1)));
  assert((Mint::UMod() - 1) % z == 0);

  if (std::min(n, m) <= 60) {
    return internal::ConvolutionNaive(lhs, rhs);
  }
  return internal::ConvolutionNTT(lhs, rhs);
}

template <typename T>
std::vector<int64_t> Convolution64Bit(const std::vector<T>& lhs,
                                      const std::vector<T>& rhs) {
  const auto n = static_cast<int>(lhs.size());
  const auto m = static_cast<int>(rhs.size());
  if (std::min(n, m) == 0) {
    return {};
  }

  // 4e18
  constexpr int64_t kMod = 4179340454199820289;
  using Mint = modular::MLong<kMod>;

  const auto z =
      static_cast<int>(std::bit_ceil(static_cast<uint32_t>(n + m - 1)));
  assert((Mint::UMod() - 1) % z == 0);

  std::vector<Mint> lhs_copy(n);
  std::copy(lhs.cbegin(), lhs.cend(), lhs_copy.begin());
  std::vector<Mint> rhs_copy(n);
  std::copy(rhs.cbegin(), rhs.cend(), rhs_copy.begin());

  const std::vector<Mint> c =
      Convolution(std::move(lhs_copy), std::move(rhs_copy));
  std::vector<int64_t> result(n + m - 1);
  for (int i = 0; i < n + m - 1; ++i) {
    result[i] = static_cast<int64_t>(c[i]);
  }
}

template <typename T>
std::vector<int> ConvolutionArbitraryMod(const std::vector<T>& lhs,
                                         const std::vector<T>& rhs, int mod) {
  const auto n = static_cast<int>(lhs.size());
  const auto m = static_cast<int>(rhs.size());
  if (std::min(n, m) == 0) {
    return {};
  }

  constexpr int64_t kMod = 4179340454199820289;
  using Mint = modular::MLong<kMod>;

  std::vector<Mint> lhs_low(n);
  std::vector<Mint> lhs_high(n);
  std::vector<Mint> lhs_sum(n);
  for (int i = 0; i < n; ++i) {
    auto [high, low] = std::div(static_cast<int>(lhs[i]), 1 << 15);
    lhs_low[i] = Mint::Raw(low);
    lhs_high[i] = Mint::Raw(high);
    lhs_sum[i] = Mint::Raw(low + high);
  }

  std::vector<Mint> rhs_low(m);
  std::vector<Mint> rhs_high(m);
  std::vector<Mint> rhs_sum(m);
  for (int i = 0; i < m; ++i) {
    auto [high, low] = std::div(static_cast<int>(rhs[i]), 1 << 15);
    rhs_low[i] = Mint::Raw(low);
    rhs_high[i] = Mint::Raw(high);
    rhs_sum[i] = Mint::Raw(low + high);
  }

  const std::vector<Mint> low_c =
      Convolution(std::move(lhs_low), std::move(rhs_low));
  const std::vector<Mint> high_c =
      Convolution(std::move(lhs_high), std::move(rhs_high));
  const std::vector<Mint> sum_c =
      Convolution(std::move(lhs_sum), std::move(rhs_sum));

  std::vector<int> result(n + m - 1);
  for (int i = 0; i < n + m - 1; ++i) {
    int64_t low = static_cast<int64_t>(low_c[i]) % mod;
    int64_t high = static_cast<int64_t>(high_c[i]) % mod;
    int64_t middle = static_cast<int64_t>(sum_c[i]) % mod;
    if (middle -= low; middle < 0) {
      middle += mod;
    }
    if (middle -= high; middle < 0) {
      middle += mod;
    }
    result[i] = static_cast<int>(low);
    if (result[i] += static_cast<int>((middle << 15) % mod); result[i] >= mod) {
      result[i] -= mod;
    }
    if (result[i] += static_cast<int>((high << 30) % mod); result[i] >= mod) {
      result[i] -= mod;
    }
  }

  return result;
}

}  // namespace ntt
