#include <bits/stdc++.h>

namespace {

template <class Fun>
class y_combinator_result {
  Fun fun_;

 public:
  template <class T>
  explicit y_combinator_result(T&& fun) : fun_(std::forward<T>(fun)) {}

  template <class... Args>
  decltype(auto) operator()(Args&&... args) {
    return fun_(std::ref(*this), std::forward<Args>(args)...);
  }
};

template <class Fun>
decltype(auto) y_combinator(Fun&& fun) {
  return y_combinator_result<std::decay_t<Fun>>(std::forward<Fun>(fun));
}

#include <vector>

// Basic cache-friendly implementation of sparse table
// $Op$ might be just a pointer to function, if you want to pass lambda, you can
// do something like:
// auto my_lambda = [](const T& lhs, const T& rhs) -> T { ... };
// SparseTable<T, static_cast<T(*)(T, T)>(my_lambda)> sparse_table(...);
// Remember that it works for idempodent binary operations only
// Time: $O(n\log n)$/$O(1)$, memory usage is $O(n\log n)
template <typename T, T (*Op)(T, T)>
class SparseTable {
 public:
  explicit SparseTable(const std::vector<T>& values) {
    const int n = static_cast<int>(values.size());  // n = 0 is not allowed
    const int log = std::__lg(n * 2);
    matrix_.resize(log);
    for (int i = 0; i < log; ++i) {
      matrix_[i].resize(n - (1 << i) + 1);
    }
    std::copy(values.cbegin(), values.cend(), matrix_.front().begin());
    for (int i = 1; i < log; ++i) {
      for (int j = 0; j < n - (1 << i) + 1; ++j) {
        matrix_[i][j] =
            Op(matrix_[i - 1][j], matrix_[i - 1][j + (1 << (i - 1))]);
      }
    }
  }

  // [first, last) : first == last is not allowed
  [[nodiscard]] T Get(int first, int last) const {
    const auto level = std::__lg(last - first);
    return Op(matrix_[level][first], matrix_[level][last - (1 << level)]);
  }

 private:
  std::vector<std::vector<T>> matrix_;
};

inline int Min(int lhs, int rhs) { return lhs < rhs ? lhs : rhs; }

// https://judge.yosupo.jp/problem/staticrmq
// https://judge.yosupo.jp/submission/179338
void RunCase([[maybe_unused]] int testcase) {
  int n;
  int q;
  std::cin >> n >> q;

  std::vector<int> v(n);
  for (int& i : v) {
    std::cin >> i;
  }
  SparseTable<int, Min> sparse_table(v);

  for (int i = 0; i < q; ++i) {
    int left;
    int right;
    std::cin >> left >> right;
    std::cout << sparse_table.Get(left, right) << "\n";
  }
}

void Main() {
  int testcases = 1;
  // std::cin >> testcases;
  for (int tt = 1; tt <= testcases; ++tt) {
    RunCase(tt);
  }
}

}  // namespace

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);
  Main();
  return 0;
}
