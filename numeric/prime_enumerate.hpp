#include <vector>

namespace numeric {

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

}  // namespace numeric
