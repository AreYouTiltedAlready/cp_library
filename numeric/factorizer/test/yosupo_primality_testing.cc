// https://judge.yosupo.jp/problem/primality_test
// https://judge.yosupo.jp/submission/180865

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

namespace math {

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

}  // namespace internal

class Factorizer {
 public:
  explicit Factorizer(int n) : n_(n), spf_(n + 1) {
    primes_.reserve(2 * n / std::__lg(n + 1));
    std::iota(spf_.begin(), spf_.end(), 0);

    for (int i = 2; i <= n; ++i) {
      if (spf_[i] == i) {
        primes_.push_back(i);
      }
      for (int p : primes_) {
        if (p > spf_[i]) {
          break;
        }
        int64_t next = static_cast<int64_t>(p) * i;
        if (next > n) {
          break;
        }
        spf_[next] = p;
      }
    }
  }

  template <typename T>
  std::map<T, int> Factorize(T n) {
    if (n == 1) {
      return {};
    }
    if (IsPrime(n)) {
      return {{n, 1}};
    }
    if (n <= n_) {
      std::map<T, int> res;
      while (n != 1) {
        res[spf_[n]] += 1;
        n /= spf_[n];
      }
      return res;
    }
    T divisor = PollardRho(n);
    auto left = Factorize(divisor);
    auto right = Factorize(n / divisor);
    if (left.size() < right.size()) {
      std::swap(left, right);
    }
    for (const auto& [p, q] : right) {
      left[p] += q;
    }
    return left;
  }

  static bool IsPrime(int n) { return IsPrime<int>(n, {2, 7, 61}); }
  static bool IsPrime(int64_t n) {
    return IsPrime<int64_t>(n,
                            {2, 325, 9375, 28178, 450775, 9780504, 1795265022});
  }

 private:
  template <typename T>
  static T Gcd(T a, T b) {
    if (a == 0 || b == 0) {
      return a + b;
    }
    int common = std::countr_zero(a | b);
    b >>= std::countr_zero(b);
    do {
      a >>= std::countr_zero(a);
      if (a < b) {
        std::swap(a, b);
      }
      a -= b;
    } while (a != 0);
    return b << common;
  }

  static inline std::mt19937_64 rng = std::mt19937_64(
      std::chrono::steady_clock::now().time_since_epoch().count());

  template <typename T>
  T PollardRho(T n) {
    for (T p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29}) {
      if (n % p == 0) {
        return p;
      }
    }

    using unsigned_t = std::make_unsigned_t<T>;
    const internal::montgomery::MontgomerySpace<T> space(n);

    unsigned_t increment{};
    auto g = [&](unsigned_t x) -> unsigned_t {
      return space.Sum(space.Product(x, x), increment);
    };

    const auto jump = static_cast<int>((sqrtl(logl(n) * sqrtl(sqrtl(n)))));
    while (true) {
      increment = space.Transform(rng() % n);
      unsigned_t start = space.Transform(rng() % n);
      unsigned_t x = start;
      unsigned_t y = start;
      unsigned_t result_gcd{};
      std::vector<unsigned_t> products(jump + 1);
      do {
        products[0] = space.Transform(1U);
        for (int i = 1; i <= jump; ++i) {
          x = g(x);
          y = g(g(y));
          products[i] = space.Product(products[i - 1], x < y ? y - x : x - y);
        }
      } while ((result_gcd = Gcd(space.Reduce(products[jump]),
                                 static_cast<unsigned_t>(n))) == 1);
      if (result_gcd == n) {
        assert(space.AreEqual(products.back(), 0));
        int index = jump;
        while (index > 0 && space.AreEqual(products[index], 0)) {
          index -= 1;
        }
        result_gcd =
            Gcd(space.Reduce(products[index]), static_cast<unsigned_t>(n));
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

    // safe because n is odd
    const internal::montgomery::MontgomerySpace<T> space(n);
    const int rank =
        std::countr_zero(static_cast<std::make_unsigned_t<T>>(n - 1));
    const T d = (n - 1) >> rank;
    const T unit = space.Transform(1);
    const T neg_unit = space.Transform(n - 1);

    for (T base : bases) {
      T base_mod = base % n;
      if (base_mod == 0) {
        continue;
      }
      base_mod = space.Power(space.Transform(base_mod), d);
      if (space.AreEqual(base_mod, unit)) {
        continue;
      }
      bool witness = true;
      for (int i = 0; i < rank; ++i) {
        if (space.AreEqual(base_mod, neg_unit)) {
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

  int n_;
  std::vector<int> spf_;
  std::vector<int> primes_;
};

}  // namespace math

void RunCase([[maybe_unused]] int testcase) {
  int64_t n;
  std::cin >> n;
  std::cout << (math::Factorizer::IsPrime(n) ? "Yes\n" : "No\n");
}

void Main() {
  int testcases = 1;
  std::cin >> testcases;
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
