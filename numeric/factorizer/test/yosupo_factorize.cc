// Problem: https://judge.yosupo.jp/problem/factorize
// Submission: https://judge.yosupo.jp/submission/181699

#include <bits/stdc++.h>
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

template <typename T>
concept unsigned_int_or_int64 =
    std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>;

namespace modular {

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

template <typename T>
concept unsigned_int_or_int64 =
    std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>;

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
        mod_inverse_(1),
        unit_(),
        neg_unit_() {
    r_square_ = static_cast<promoted_t>(r_square_) * r_square_ % mod;
    for (int i = 0; i < 6; ++i) {
      mod_inverse_ *= static_cast<unsigned_t>(2) - mod_ * mod_inverse_;
    }
    unit_ = Transform(1);
    neg_unit_ = Transform(mod_ - 1);
  }

  [[nodiscard]] constexpr unsigned_t Unit() const { return unit_; }

  [[nodiscard]] constexpr unsigned_t NegUnit() const { return neg_unit_; }

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

  [[nodiscard]] constexpr bool TestEquality(unsigned_t lhs,
                                            unsigned_t rhs) const {
    return (lhs < mod_ ? lhs : lhs - mod_) == (rhs < mod_ ? rhs : rhs - mod_);
  }

  [[nodiscard]] constexpr unsigned_t Power(unsigned_t x, uint64_t n) const {
    unsigned_t result = unit_;
    while (n > 0) {
      if (n % 2 == 1) {
        result = Product(result, x);
      }
      x = Product(x, x);
      n /= 2;
    }
    return result;
  }

  [[nodiscard]] constexpr unsigned_t ExactReduce(promoted_t n) const {
    unsigned_t result = Reduce(n);
    return result - (result >= mod_ ? mod_ : 0);
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
  unsigned_t unit_;
  unsigned_t neg_unit_;
};

}  // namespace montgomery

}  // namespace modular

namespace math {

class Factorizer {
 public:
  constexpr Factorizer() : spf_(), primes_(), magic_() {
    spf_[1] = 1;
    for (int i = 2; i <= kSieveBound; ++i) {
      if (spf_[i] == 0) {
        spf_[i] = i;
        primes_[primes_count_] = i;
        magic_[primes_count_] = static_cast<uint64_t>(-1) / i + 1;
        primes_count_ += 1;
      }
      for (int j = 0; j < primes_count_; ++j) {
        if (primes_[j] > spf_[i]) {
          break;
        }
        if (auto next = static_cast<int64_t>(primes_[j]) * i;
            next <= kSieveBound) {
          spf_[next] = primes_[j];
        } else {
          break;
        }
      }
    }
  }

  template <modular::unsigned_int_or_int64 T>
  std::vector<std::pair<T, int>> Factorize(T n) const {
    if (n == 1) {
      return {};
    }

    if (n < (1 << 30)) {
      if (IsPrime(static_cast<int>(n))) {
        return {{n, 1}};
      }

      if (n <= kSieveBound) {
        return FactorizeSmall(n);
      }
      return FactorizeMedium(n);
    }

    if (IsPrime(static_cast<int64_t>(n))) {
      return {{n, 1}};
    }

    T divisor = PollardRho(n);
    return MergeFactors(Factorize(divisor), Factorize(n / divisor));
  }

  static bool IsPrime(int n) {
    return n > 0 && IsPrime<uint32_t>(n, {2, 7, 61});
  }

  static bool IsPrime(int64_t n) {
    return n > 0 && IsPrime<uint64_t>(
                        n, {2, 325, 9375, 28178, 450775, 9780504, 1795265022});
  }

 private:
  template <std::integral T>
  [[nodiscard]] std::vector<std::pair<T, int>> FactorizeSmall(T n) const {
    std::vector<std::pair<T, int>> result{};
    while (n != 1) {
      int d = spf_[n];
      int count = 0;
      while (n % d == 0) {
        n /= d;
        count += 1;
      }
      result.emplace_back(d, count);
    }
    return result;
  }

  template <std::integral T>
  [[nodiscard]] std::vector<std::pair<T, int>> FactorizeMedium(T n) const {
    auto unsigned_n = static_cast<uint32_t>(n);
    std::vector<std::pair<T, int>> result{};
    for (int i = 0; i < primes_count_ && primes_[i] <= unsigned_n; ++i) {
      if (magic_[i] * unsigned_n >= magic_[i]) {
        continue;
      }
      int d = primes_[i];
      int count = 0;
      while (magic_[i] * unsigned_n < magic_[i]) {
        unsigned_n /= d;
        count += 1;
      }
      result.emplace_back(d, count);
    }
    if (unsigned_n != 1) {
      result.emplace_back(static_cast<T>(unsigned_n), 1);
    }
    return result;
  }

  template <std::integral T>
  static std::vector<std::pair<T, int>> MergeFactors(
      const std::vector<std::pair<T, int>>& lhs,
      const std::vector<std::pair<T, int>>& rhs) {
    auto lhs_it = lhs.cbegin();
    auto rhs_it = rhs.cbegin();
    std::vector<std::pair<T, int>> result{};
    result.reserve(lhs.size() + rhs.size());
    while (lhs_it != lhs.cend() && rhs_it != rhs.cend()) {
      if ((*lhs_it).first == (*rhs_it).first) {
        result.emplace_back((*lhs_it).first,
                            (*lhs_it).second + (*rhs_it).second);
        ++lhs_it;
        ++rhs_it;
      } else if ((*lhs_it).first < (*rhs_it).first) {
        result.push_back(*lhs_it);
        ++lhs_it;
      } else {
        result.push_back(*rhs_it);
        ++rhs_it;
      }
    }
    while (lhs_it != lhs.cend()) {
      result.push_back(*lhs_it);
      ++lhs_it;
    }
    while (rhs_it != rhs.cend()) {
      result.push_back(*rhs_it);
      ++rhs_it;
    }
    return result;
  }

  static inline std::mt19937_64 gen = std::mt19937_64(
      std::chrono::steady_clock::now().time_since_epoch().count());

  static uint64_t PollardRho(uint64_t n) {
    for (int p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29}) {
      if (n % p == 0) {
        return p;
      }
    }

    uint64_t increment{};
    const modular::montgomery::MontgomerySpace<uint64_t> space(n);
    auto g = [&](uint64_t x) -> uint64_t {
      return space.Sum(space.Product(x, x), increment);
    };

    const auto jump = static_cast<int>((sqrtl(logl(n) * sqrtl(sqrtl(n)))));
    while (true) {
      increment = space.Transform(gen() % n);
      uint64_t start = space.Transform(gen() % n);
      uint64_t x = start;
      uint64_t y = start;
      std::vector<uint64_t> products(jump + 1);

      uint64_t result_gcd = 1;
      while (result_gcd == 1) {
        products[0] = space.Unit();
        for (int i = 1; i <= jump; ++i) {
          x = g(x);
          y = g(g(y));
          products[i] = space.Product(products[i - 1], x < y ? y - x : x - y);
        }
        result_gcd = std::gcd(space.Reduce(products.back()), n);
      }

      if (result_gcd == n) {
        int index = jump;
        while (index > 0 && space.TestEquality(products[index], 0)) {
          index -= 1;
        }
        result_gcd = std::gcd(space.Reduce(products[index]), n);
      }

      if (result_gcd != n && result_gcd != 1) {
        return result_gcd;
      }
    }
  }

  template <typename T>
  static bool IsPrime(T n, const std::vector<T>& bases) {
    if (n < 2) {
      return false;
    }
    static constexpr std::array<T, 10> kSmallPrimes = {2,  3,  5,  7,  11,
                                                       13, 17, 19, 23, 29};
    for (T p : kSmallPrimes) {
      if (n % p == 0) {
        return n == p;
      }
    }

    if (n < 31 * 31) {
      return true;
    }

    using unsigned_t = std::make_unsigned_t<T>;
    const modular::montgomery::MontgomerySpace<unsigned_t> space(n);
    const int rank =
        std::countr_zero(static_cast<std::make_unsigned_t<T>>(n - 1));
    const T d = (n - 1) >> rank;

    for (T base : bases) {
      T base_mod = base % n;
      if (base_mod == 0) {
        continue;
      }
      base_mod = space.Power(space.Transform(base_mod), d);
      if (space.TestEquality(base_mod, space.Unit())) {
        continue;
      }
      bool witness = true;
      for (int i = 0; i < rank; ++i) {
        if (space.TestEquality(base_mod, space.NegUnit())) {
          witness = false;
          break;
        }
        base_mod = space.Product(base_mod, base_mod);
      }
      if (witness) {
        return false;
      }
    }

    return true;
  }

  static constexpr int kSieveBound = 1 << 16;

  std::array<int, kSieveBound + 1> spf_;
  std::array<int, kSieveBound + 1> primes_;
  std::array<uint64_t, kSieveBound + 1> magic_;
  int primes_count_{};
};

}  // namespace math

void RunCase([[maybe_unused]] int testcase) {
  math::Factorizer factorizer;
  int q{};
  std::cin >> q;
  while (q--) {
    uint64_t x{};
    std::cin >> x;
    auto result = factorizer.Factorize(x);
    int q_sum = 0;
    for (const auto& [p, q] : result) {
      q_sum += q;
    }
    std::cout << q_sum << " ";
    for (const auto& [p, q] : result) {
      for (int i = 0; i < q; ++i) {
        std::cout << p << " ";
      }
    }
    std::cout << "\n";
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
