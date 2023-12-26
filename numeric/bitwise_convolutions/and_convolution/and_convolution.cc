#include <type_traits>
#include <vector>

template <typename T>
void ZetaTransform(std::vector<T>& v) {
  const int n = static_cast<int>(v.size());
  for (int i = 1; i < n; i *= 2) {
    for (int j = 0; j < n; ++j) {
      if ((j & i) != 0) {
        v[j ^ i] += v[j];
      }
    }
  }
}

template <typename T>
void MobiusTransform(std::vector<T>& v) {
  const int n = static_cast<int>(v.size());
  for (int i = n / 2; i > 0; i /= 2) {
    for (int j = n - 1; j >= 0; --j) {
      if ((j & i) != 0) {
        v[j ^ i] -= v[j];
      }
    }
  }
}

template <typename T, typename U>
std::vector<std::common_type_t<T, U>> AndConvolution(std::vector<T> lhs,
                                                     std::vector<U> rhs) {
  const int n = static_cast<int>(lhs.size());
  ZetaTransform(lhs);
  ZetaTransform(rhs);
  using V = std::common_type_t<T, U>;
  std::vector<V> result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = static_cast<V>(lhs[i]) * rhs[i];
  }
  MobiusTransform(result);
  return result;
}
