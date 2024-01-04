#include <vector>

std::vector<int> PrimeEnumerate(int n) {
  std::vector<int> primes;
  primes.reserve(1.1 * n / std::__lg(n * 2));
  std::vector<bool> is_prime(n + 1, true);
  for (int i = 2; i <= n; ++i) {
    if (!is_prime[i]) {
      continue;
    }
    primes.push_back(i);
    for (int j = i * 2; j <= n; j += i) {
      is_prime[j] = false;
    }
  }
  return primes;
}

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
