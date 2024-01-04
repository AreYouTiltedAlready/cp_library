// Problem: https://judge.yosupo.jp/problem/primality_test
// Submission: https://judge.yosupo.jp/submission/181698

#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <type_traits>
#include <utility>
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

namespace barrett {

using uint128_t = __uint128_t;

constexpr uint128_t MultiplyHigh(uint128_t lhs, uint128_t rhs) {
  const auto a = static_cast<uint64_t>(lhs >> 64);
  const auto b = static_cast<uint64_t>(lhs);
  const auto c = static_cast<uint64_t>(rhs >> 64);
  const auto d = static_cast<uint64_t>(rhs);
  uint128_t ac = static_cast<uint128_t>(a) * c;
  uint128_t ad = static_cast<uint128_t>(a) * d;
  uint128_t bc = static_cast<uint128_t>(b) * c;
  uint128_t bd = static_cast<uint128_t>(b) * d;
  uint128_t carry = static_cast<uint128_t>(static_cast<uint64_t>(ad)) +
                    static_cast<uint128_t>(static_cast<uint64_t>(bc)) +
                    (bd >> 64);
  return ac + (ad >> 64) + (bc >> 64) + (carry >> 64);
}

template <unsigned_int_or_int64 T>
class Barrett {
 public:
  using promoted_t = integral_promotion_t<T>;

  constexpr explicit Barrett(T mod)
      : mod_inverse_(static_cast<promoted_t>(-1) / mod + 1), mod_(mod) {}

  [[nodiscard]] constexpr T Product(T lhs, T rhs) const {
    return Reduce(static_cast<promoted_t>(lhs) * rhs);
  }

  [[nodiscard]] constexpr T Reduce(promoted_t n) const {
    promoted_t m = mod_;
    if constexpr (std::is_same_v<T, uint32_t>) {
      m *= static_cast<promoted_t>((static_cast<uint128_t>(n) * mod_inverse_) >>
                                   64);
    } else {
      m *= MultiplyHigh(n, mod_inverse_);
    }
    n += (n < m ? mod_ : 0);
    return n - m;
  }

 private:
  promoted_t mod_inverse_;
  T mod_;
};

}  // namespace barrett

}  // namespace modular

namespace math {

class Factorizer {
 public:
  explicit Factorizer(int n) : spf_(n + 1), primes_(), n_(n) {
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
    const modular::barrett::Barrett<unsigned_t> barrett(n);

    unsigned_t increment{};
    auto g = [&](unsigned_t x) -> unsigned_t {
      x = barrett.Product(x, x);
      if ((x += increment) >= n) {
        x -= n;
      }
      return x;
    };

    const auto jump = static_cast<int>((sqrtl(logl(n) * sqrtl(sqrtl(n)))));
    while (true) {
      increment = rng() % n;
      unsigned_t start = rng() % n;
      unsigned_t x = start;
      unsigned_t y = start;
      unsigned_t result_gcd{};
      std::vector<unsigned_t> products(jump + 1);
      do {
        products[0] = 1;
        for (int i = 1; i <= jump; ++i) {
          x = g(x);
          y = g(g(y));
          products[i] = barrett.Product(products[i - 1], x < y ? y - x : x - y);
        }
      } while ((result_gcd =
                    std::gcd(products[jump], static_cast<unsigned_t>(n))) == 1);
      if (result_gcd == n) {
        assert(products.back() == 0);
        int index = jump;
        while (index > 0 && products[index] == 0) {
          index -= 1;
        }
        result_gcd = std::gcd(products[index], static_cast<unsigned_t>(n));
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
    const modular::barrett::Barrett<unsigned_t> barrett(n);
    const int rank =
        std::countr_zero(static_cast<std::make_unsigned_t<T>>(n - 1));
    const T d = (n - 1) >> rank;

    auto Power = [&barrett](unsigned_t x, unsigned_t n) -> unsigned_t {
      unsigned_t res = 1;
      while (n > 0) {
        if (n % 2 == 1) {
          res = barrett.Product(res, x);
        }
        x = barrett.Product(x, x);
        n /= 2;
      }
      return res;
    };

    for (T base : bases) {
      T base_mod = base % n;
      if (base_mod == 0) {
        continue;
      }
      base_mod = Power(base_mod, d);
      if (base_mod == 1) {
        continue;
      }
      bool witness = true;
      for (int i = 0; i < rank; ++i) {
        if (base_mod == n - 1) {
          witness = false;
          break;
        }
        base_mod = barrett.Product(base_mod, base_mod);
      }
      if (witness) {
        return false;
      }
    }

    return true;
  }

  std::vector<int> spf_;
  std::vector<int> primes_;

  int n_;
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
