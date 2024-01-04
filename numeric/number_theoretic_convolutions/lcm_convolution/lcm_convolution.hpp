// #include "numeric/prime_enumerate.hpp"

#include <type_traits>
#include <vector>

namespace numeric {
namespace lcm_convolution {

template <typename T>
void ZetaTransform(std::vector<T>& v, const std::vector<int>& primes) {
  const int n = static_cast<int>(v.size()) - 1;
  for (int p : primes) {
    for (int i = 1; i <= n / p; ++i) {
      v[i * p] += v[i];
    }
  }
}

template <typename T>
void MobiusTransform(std::vector<T>& v, const std::vector<int>& primes) {
  const int n = static_cast<int>(v.size()) - 1;
  for (int p : primes) {
    for (int i = n / p; i > 0; --i) {
      v[i * p] -= v[i];
    }
  }
}

template <typename T, typename U>
std::vector<std::common_type_t<T, U>> LcmConvolution(std::vector<T> lhs,
                                                     std::vector<U> rhs) {
  const int n = static_cast<int>(lhs.size()) - 1;
  const std::vector<int> primes = PrimeEnumerate(n);

  ZetaTransform(lhs, primes);
  ZetaTransform(rhs, primes);
  using V = std::common_type_t<T, U>;
  std::vector<V> result(n + 1);
  for (int i = 0; i <= n; ++i) {
    result[i] = static_cast<V>(lhs[i]) * rhs[i];
  }
  MobiusTransform(result, primes);
  return result;
}

}  // namespace lcm_convolution
}  // namespace numeric
