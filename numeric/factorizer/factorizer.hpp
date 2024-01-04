// #include "modular/barrett.hpp"

#include <array>
#include <bit>
#include <chrono>
#include <map>
#include <numeric>
#include <random>
#include <vector>

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
