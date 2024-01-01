// https://judge.yosupo.jp/problem/convolution_mod
// https://judge.yosupo.jp/submission/180867

#include <array>
#include <bit>
#include <cassert>
#include <chrono>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <vector>

namespace {

template <class Fun>
class y_combinator_result {
  Fun fun_;

 public:
  template <class T>
  explicit y_combinator_result(T&& fun)  // NOLINT(*forwarding-reference*)
      : fun_(std::forward<T>(fun)) {}

  template <class... Args>
  decltype(auto) operator()(Args&&... args) {
    return fun_(std::ref(*this), std::forward<Args>(args)...);
  }
};

template <class Fun>
decltype(auto) y_combinator(Fun&& fun) {
  return y_combinator_result<std::decay_t<Fun>>(std::forward<Fun>(fun));
}

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

namespace ntt {

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
  int d = (mod - 1) >> std::countr_zero(static_cast<uint32_t>(mod - 1));

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

    using Mint = modular::MInt<kMod>;

    const int z =
        static_cast<int>(std::bit_ceil(static_cast<uint32_t>((n + m - 1))));
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

  template <modular::internal::ModInt Mint>
  std::vector<Mint> Convolution(const std::vector<Mint>& lhs,
                                const std::vector<Mint>& rhs) {
    const int n = static_cast<int>(lhs.size());
    const int m = static_cast<int>(rhs.size());
    if (std::min(n, m) == 0) {
      return {};
    }

    const int z =
        static_cast<int>(std::bit_ceil(static_cast<uint32_t>((n + m - 1))));
    assert((Mint::UMod() - 1) % z == 0);

    if (std::min(n, m) <= 60) {
      return ConvolutionNaive(lhs, rhs);
    }

    return ConvolutionNTT(lhs, rhs);
  }

  template <modular::internal::ModInt Mint>
  std::vector<Mint> Convolution(std::vector<Mint>&& lhs,
                                std::vector<Mint>&& rhs) {
    const int n = static_cast<int>(lhs.size());
    const int m = static_cast<int>(rhs.size());
    if (std::min(n, m) == 0) {
      return {};
    }

    const int z =
        static_cast<int>(std::bit_ceil(static_cast<uint32_t>((n + m - 1))));
    assert((Mint::UMod() - 1) % z == 0);

    if (std::min(n, m) <= 60) {
      return ConvolutionNaive(lhs, rhs);
    }

    return ConvolutionNTT(lhs, rhs);
  }

 private:
  template <modular::internal::ModInt Mint,
            int kG = kPrimitiveRoot<Mint::UMod()>>
  struct NTTInfo {
    static constexpr int kBinaryRank = std::countr_zero(Mint::UMod() - 1);

    constexpr NTTInfo() : roots() {
      roots[kBinaryRank - 1] =
          Power(Mint(kG), (Mint::UMod() - 1) >> kBinaryRank);
      for (int i = kBinaryRank - 2; i >= 0; --i) {
        roots[i] = roots[i + 1] * roots[i + 1];
      }
    }

    std::array<Mint, kBinaryRank> roots;
  };

  template <modular::internal::ModInt Mint>
  std::vector<Mint> ConvolutionNTT(std::vector<Mint> lhs,
                                   std::vector<Mint> rhs) {
    const int n = static_cast<int>(lhs.size());
    const int m = static_cast<int>(rhs.size());
    const int z =
        static_cast<int>(std::bit_ceil(static_cast<uint32_t>((n + m - 1))));
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

  template <modular::internal::ModInt Mint>
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

  template <modular::internal::ModInt Mint>
  void Butterfly(std::vector<Mint>& data) {
    const int n = static_cast<int>(data.size());
    const int rank = std::countr_zero(static_cast<uint32_t>(n));

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

  std::array<int, 1 << kMaxRank> bit_inverse_;
};

}  // namespace ntt

void RunCase([[maybe_unused]] int testcase) {
  int n;
  int m;
  std::cin >> n >> m;

  using Mint = modular::MInt<998244353>;
  std::vector<Mint> p(n);
  for (Mint& it : p) {
    std::cin >> it;
  }

  std::vector<Mint> q(m);
  for (Mint& it : q) {
    std::cin >> it;
  }

  ntt::NTT<20> ntt{};
  std::vector<Mint> r = ntt.Convolution(p, q);
  for (const Mint& it : r) {
    std::cout << it << " ";
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
