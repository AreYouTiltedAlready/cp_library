// #include "modular/modint/static_modint.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>

namespace numeric {
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
  int d = static_cast<int>((mod - 1) >>
                           std::countr_zero(static_cast<uint32_t>(mod - 1)));

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

template <modular::modint::modint Mint,
          int kG = kPrimitiveRoot<static_cast<int>(Mint::UMod())>>
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

template <modular::modint::modint Mint>
void Butterfly(std::vector<Mint>& data) {
  const auto n = static_cast<int>(data.size());
  const int rank = std::countr_zero(static_cast<uint32_t>(n));

  static constexpr NTTInfo<Mint> info{};
  using unsigned_t = modular::integral_promotion_t<typename Mint::unsigned_t>;

  int length = 0;  // data[i, i + (n >> length), i + 2 * (n >> length), ..] is
                   // transformed
  while (length < rank) {
    if (rank - length == 1) {
      const int p = 1 << (rank - length - 1);
      Mint root = Mint::Raw(1U);
      for (int s = 0; s < (1 << length); s++) {
        const int offset = s << (rank - length);
        for (int i = 0; i < p; i++) {
          Mint left = data[i + offset];
          Mint right = data[i + offset + p] * root;
          data[i + offset] = left + right;
          data[i + offset + p] = left - right;
        }
        if (s + 1 != (1 << length)) {
          root *= info.rate_2[std::countr_one(static_cast<uint32_t>(s))];
        }
      }
      length += 1;
    } else {
      const int p = 1 << (rank - length - 2);
      Mint root = Mint::Raw(1U);
      const Mint imag = info.root[2];
      for (int s = 0; s < (1 << length); s++) {
        const Mint root_square = root * root;
        const Mint root_cube = root * root_square;
        const int offset = s << (rank - length);
        for (int i = 0; i < p; i++) {
          const unsigned_t squared_mod =
              static_cast<unsigned_t>(Mint::UMod()) * Mint::UMod();
          const auto a_0 = static_cast<unsigned_t>(data[i + offset]);
          const auto a_1 =
              static_cast<unsigned_t>(data[i + offset + p]) * root.Get();
          const auto a_2 = static_cast<unsigned_t>(data[i + offset + 2 * p]) *
                           root_square.Get();
          const auto a_3 = static_cast<unsigned_t>(data[i + offset + 3 * p]) *
                           root_cube.Get();
          const auto aux =
              static_cast<unsigned_t>(Mint(a_1 + squared_mod - a_3)) *
              imag.Get();
          const auto neg_a_2 = squared_mod - a_2;
          data[i + offset] = a_0 + a_2 + a_1 + a_3;
          data[i + offset + 1 * p] = a_0 + a_2 + (squared_mod * 2 - a_1 - a_3);
          data[i + offset + 2 * p] = a_0 + neg_a_2 + aux;
          data[i + offset + 3 * p] = a_0 + neg_a_2 + squared_mod - aux;
        }
        if (s + 1 != (1 << length)) {
          root *= info.rate_3[std::countr_one(static_cast<uint32_t>(s))];
        }
      }
      length += 2;
    }
  }
}

template <modular::modint::modint Mint>
void InverseButterfly(std::vector<Mint>& data) {
  const auto n = static_cast<int>(data.size());
  const int rank = std::countr_zero(static_cast<uint32_t>(n));

  static constexpr NTTInfo<Mint> info{};
  using unsigned_t = modular::integral_promotion_t<typename Mint::unsigned_t>;

  int length = rank;  // data[i, i + (n >> length), i + 2 * (n >> length, ..] is
                      // transformed
  while (length > 0) {
    if (length == 1) {
      const int p = 1 << (rank - length);
      Mint inv_root = Mint::Raw(1U);
      for (int s = 0; s < (1 << (length - 1)); s++) {
        const int offset = s << (rank - length + 1);
        for (int i = 0; i < p; i++) {
          Mint left = data[i + offset];
          Mint right = data[i + offset + p];
          data[i + offset] = left + right;
          data[i + offset + p] = (left - right) * inv_root;
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
          const auto a_0 = static_cast<unsigned_t>(data[i + offset]);
          const auto a_1 = static_cast<unsigned_t>(data[i + offset + 1 * p]);
          const auto a_2 = static_cast<unsigned_t>(data[i + offset + 2 * p]);
          const auto a_3 = static_cast<unsigned_t>(data[i + offset + 3 * p]);
          const auto aux = static_cast<unsigned_t>(
              Mint((a_2 + Mint::UMod() - a_3) * inv_imag.Get()));
          data[i + offset] = a_0 + a_1 + a_2 + a_3;
          data[i + offset + 1 * p] =
              (a_0 + Mint::UMod() - a_1 + aux) * inv_root.Get();
          data[i + offset + 2 * p] =
              (a_0 + a_1 + Mint::UMod() - a_2 + Mint::UMod() - a_3) *
              inv_root_square.Get();
          data[i + offset + 3 * p] =
              (a_0 + Mint::UMod() - a_1 + Mint::UMod() - aux) *
              inv_root_cube.Get();
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

template <modular::modint::modint Mint>
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

template <modular::modint::modint Mint>
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

template <typename T>
constexpr std::array<T, 3> ExtendedGcd(T a, T b) {
  std::array x{T{1}, T{0}};
  while (b > 0) {
    x[0] = std::exchange(x[1], x[0] - a / b * x[1]);
    a = std::exchange(b, a % b);
  }
  return {a, x[0], x[1]};
}

template <modular::modint::modint Mint>
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

  if (std::min(n, m) <= 60) {
    return internal::ConvolutionNaive(lhs, rhs);
  }
  return internal::ConvolutionNTT(std::move(lhs), std::move(rhs));
}

template <modular::modint::modint Mint>
std::vector<Mint> Convolution(const std::vector<Mint>& lhs,
                              const std::vector<Mint>& rhs) {
  if (lhs.empty() || rhs.empty()) {
    return {};
  }
  return Convolution(std::move(lhs), std::move(rhs));
}

template <int kMod = 998244353, typename T>
std::vector<modular::modint::StaticMInt<kMod>> ConvolutionMod(
    const std::vector<T>& lhs, const std::vector<T>& rhs) {
  using Mint = modular::modint::StaticMInt<kMod>;

  const auto n = static_cast<int>(lhs.size());
  const auto m = static_cast<int>(rhs.size());
  if (std::min(n, m) == 0) {
    return {};
  }
  const auto z =
      static_cast<int>(std::bit_ceil(static_cast<uint32_t>(n + m - 1)));
  assert((Mint::UMod() - 1) % z == 0);

  std::vector<Mint> lhs_copy(n);
  for (int i = 0; i < n; ++i) {
    lhs_copy[i] = static_cast<int64_t>(lhs[i]);
  }

  std::vector<Mint> rhs_copy(m);
  for (int i = 0; i < m; ++i) {
    rhs_copy[i] = static_cast<int64_t>(rhs[i]);
  }

  return Convolution(std::move(lhs_copy), std::move(rhs_copy));
}

struct MultiModConvolutionResult {
  static constexpr std::array<int, 3> kMods = {
      {754974721, 167772161, 469762049}};

  static constexpr int64_t kProduct_01 =
      static_cast<int64_t>(kMods[0]) * kMods[1];

  static constexpr int kInverse_01 = static_cast<int>(
      modular::modint::StaticMInt<kMods[1]>(kMods[0]).Inverse());
  static constexpr int kInverse_10 = static_cast<int>(
      modular::modint::StaticMInt<kMods[0]>(kMods[1]).Inverse());
  static constexpr int kInverse_01_2 =
      static_cast<int>(modular::modint::StaticMInt<kMods[2]>(
                           static_cast<int64_t>(kMods[0]) * kMods[1])
                           .Inverse());

  using Mint_0 = modular::modint::StaticMInt<kMods[0]>;
  using Mint_1 = modular::modint::StaticMInt<kMods[1]>;
  using Mint_2 = modular::modint::StaticMInt<kMods[2]>;

  explicit MultiModConvolutionResult(std::vector<Mint_0>&& c_0,
                                     std::vector<Mint_1>&& c_1,
                                     std::vector<Mint_2>&& c_2)
      : c_0(std::move(c_0)), c_1(std::move(c_1)), c_2(std::move(c_2)) {}

  std::vector<int> Restore(int mod) {
    std::vector<int> result(c_0.size());
    for (int i = 0; i < static_cast<int>(c_0.size()); ++i) {
      auto aux_10 = static_cast<int64_t>(
          (static_cast<__int128_t>(kMods[1]) * c_0[i].Get() % kProduct_01) *
          kInverse_10 % kProduct_01);
      auto aux_01 = static_cast<int64_t>(
          (static_cast<__int128_t>(kMods[0]) * c_1[i].Get() % kProduct_01) *
          kInverse_01 % kProduct_01);
      int64_t aux_sum = aux_01 + aux_10;
      if (aux_sum >= kProduct_01) {
        aux_sum -= kProduct_01;
      }
      int64_t c = aux_sum;
      aux_sum = static_cast<int64_t>(c_2[i]) - (aux_sum % kMods[2]);
      if (aux_sum < 0) {
        aux_sum += kMods[2];
      }
      aux_sum *= kInverse_01_2;
      aux_sum %= kMods[2];
      int64_t value = ((kProduct_01 % mod) * aux_sum % mod);
      if ((value += c % mod) >= mod) {
        value -= mod;
      }
      result[i] = static_cast<int>(value);
    }
    return result;
  }

  std::vector<Mint_0> c_0;
  std::vector<Mint_1> c_1;
  std::vector<Mint_2> c_2;
};

template <typename T>
MultiModConvolutionResult MultiModConvolution(const std::vector<T>& lhs,
                                              const std::vector<T>& rhs) {
  const auto n = static_cast<int>(lhs.size());
  const auto m = static_cast<int>(rhs.size());
  if (std::min(n, m) == 0) {
    return MultiModConvolutionResult({}, {}, {});
  }

  const auto z =
      static_cast<int>(std::bit_ceil(static_cast<uint32_t>(n + m - 1)));
  assert(z <= (1 << 24));

  return MultiModConvolutionResult(
      ConvolutionMod<MultiModConvolutionResult::kMods[0]>(lhs, rhs),
      ConvolutionMod<MultiModConvolutionResult::kMods[1]>(lhs, rhs),
      ConvolutionMod<MultiModConvolutionResult::kMods[2]>(lhs, rhs));
}

template <typename T>
std::vector<int> ConvolutionArbitraryMod(const std::vector<T>& lhs,
                                         const std::vector<T>& rhs, int mod) {
  auto result = MultiModConvolution(lhs, rhs);
  return result.Restore(mod);
}

}  // namespace ntt
}  // namespace numeric
