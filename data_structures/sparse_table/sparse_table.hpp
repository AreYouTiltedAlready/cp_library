#include <cstdint>
#include <vector>

namespace ds {

namespace sparse_table {
// Basic cache-friendly implementation of sparse table
// $Op$ might be just a pointer to function, if you want to pass lambda, you can
// do something like:
// auto my_lambda = [](T lhs, T rhs) -> T { ... };
// SparseTable<T, static_cast<T(*)(T, T)>(my_lambda)> sparse_table(...);
// Remember that it works for idempotent binary operations only
// Time: $O(n\log n)$/$O(1)$, memory usage is $O(n\log n)
template <typename T, auto Op>
class SparseTable {
  static_assert(
      std::is_convertible_v<decltype(Op), T (*)(T, T)> ||
          std::is_convertible_v<decltype(Op), T (*)(const T&, const T&)>,
      "Op must work as S(S, S)");

 public:
  template <typename U,
            std::enable_if_t<std::is_same_v<std::decay_t<U>, std::vector<T>>,
                             void>* = nullptr>
  explicit SparseTable(U&& values) {
    const int n = static_cast<int>(values.size());  // n = 0 is not allowed
    const int log = std::bit_width(static_cast<uint32_t>(n));
    matrix_.resize(log);
    for (int i = 0; i < log; ++i) {
      matrix_[i].resize(n - (1 << i) + 1);
    }
    matrix_[0] = std::forward<U>(values);
    for (int i = 1; i < log; ++i) {
      for (int j = 0; j < n - (1 << i) + 1; ++j) {
        matrix_[i][j] =
            Op(matrix_[i - 1][j], matrix_[i - 1][j + (1 << (i - 1))]);
      }
    }
  }

  // [first, last) : first == last is not allowed
  [[nodiscard]] T Get(int first, int last) const {
    const int level = std::bit_width(static_cast<uint32_t>(last - first)) - 1;
    return Op(matrix_[level][first], matrix_[level][last - (1 << level)]);
  }

 private:
  std::vector<std::vector<T>> matrix_;
};

}  // namespace sparse_table
}  // namespace ds
