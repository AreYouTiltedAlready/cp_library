// #include "numeric/prime_enumerate.hpp"

#include <vector>

namespace numeric {
namespace gcd_convolution {

template <typename T>
void ZetaTransform(std::vector<T>& v, const std::vector<int>& primes) {
  const int n = static_cast<int>(v.size()) - 1;
  for (int p : primes) {
    for (int i = n / p; i > 0; --i) {
      v[i] += v[i * p];
    }
  }
}

template <typename T>
void MobiusTransform(std::vector<T>& v, const std::vector<int>& primes) {
  const int n = static_cast<int>(v.size()) - 1;
  for (int p : primes) {
    for (int i = 1; i <= n / p; ++i) {
      v[i] -= v[i * p];
    }
  }
}

template <typename T, typename U>
std::vector<std::common_type_t<T, U>> GcdConvolution(std::vector<T> lhs,
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

}  // namespace gcd_convolution
}  // namespace numeric
