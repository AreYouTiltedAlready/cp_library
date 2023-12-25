#include <type_traits>
#include <vector>

template <typename T>
void FWHT(std::vector<T>& v) {
  const int n = static_cast<int>(v.size());
  for (int i = 1; i < n; i *= 2) {
    for (int j = 0; j < n; j += i * 2) {
      for (int k = 0; k < i; ++k) {
        T x = v[j + k];
        T y = v[i + j + k];
        v[j + k] = x + y;
        v[i + j + k] = x - y;
      }
    }
  }
}

template <typename T, typename U>
std::vector<std::common_type_t<T, U>> XorConvolution(std::vector<T> lhs,
                                                     std::vector<U> rhs) {
  const int n = static_cast<int>(lhs.size());
  FWHT(lhs);
  FWHT(rhs);
  using V = std::common_type_t<T, U>;
  std::vector<V> result(n);
  for (int i = 0; i < n; ++i) { result[i] = static_cast<V>(lhs[i]) * rhs[i]; }
  FWHT(result);
  for (V& c : result) { c /= n; }
  return result;
}