#include <numeric>
#include <type_traits>
#include <vector>

// Just a version of sparse table for Range M(ax)|(in)imum Query problem
// GetIndex(first, last) always yields the FIRST occurrence of min/max in range
// Note that is has not default constructor (it's quite problematic to provide
// it)
// So, for stuff like fast lca class, which has to hold an instance of rmq
// solver, use pointer/optional/whatever..

enum class RmqMode {
  kMax,
  kMin,
};

template <typename T, RmqMode mode>
class RmqSolver {
 public:
  template <typename U, typename std::enable_if_t<
                            std::is_same_v<std::decay_t<U>, std::vector<T>>,
                            void>* = nullptr>
  explicit RmqSolver(U values) : values_(values) {
    const int n = static_cast<int>(values.size());  // n = 0 is not allowed
    const int log = std::__lg(n * 2);
    matrix_.resize(log);
    for (int i = 0; i < log; ++i) {
      matrix_[i].resize(n - (1 << i) + 1);
    }
    std::iota(matrix_.front().begin(), matrix_.front().end(), 0);
    for (int i = 1; i < log; ++i) {
      for (int j = 0; j < n - (1 << i) + 1; ++j) {
        matrix_[i][j] =
            Merger(matrix_[i - 1][j], matrix_[i - 1][j + (1 << (i - 1))]);
      }
    }
  }

  // [first, last) : first == last is not allowed
  [[nodiscard]] T GetIndex(int first, int last) const {
    const auto level = std::__lg(last - first);
    return Merger(matrix_[level][first], matrix_[level][last - (1 << level)]);
  }

  [[nodiscard]] T GetValue(int first, int last) const {
    return values_[GetIndex(first, last)];
  }

 private:
  inline int Merger(int lhs, int rhs) const {
    if constexpr (mode == RmqMode::kMin) {
      return values_[lhs] < values_[rhs] ? lhs : rhs;
    }
    return values_[lhs] > values_[rhs] ? lhs : rhs;
  }

  std::vector<T> values_;
  std::vector<std::vector<int>> matrix_;
};

template <typename T>
using RmqSolverMax = RmqSolver<T, RmqMode::kMax>;

template <typename T>
using RmqSolverMin = RmqSolver<T, RmqMode::kMin>;
