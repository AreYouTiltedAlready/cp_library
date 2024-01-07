// #include "modular/montgomery.hpp"

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