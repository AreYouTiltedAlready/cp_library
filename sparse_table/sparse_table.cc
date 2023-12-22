#include <iostream>
#include <vector>

template <typename T, auto Op>
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

// https://judge.yosupo.jp/problem/staticrmq
// https://judge.yosupo.jp/submission/178572

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);

  int n;
  int q;
  std::cin >> n >> q;

  std::vector<int> v(n);
  for (int& i : v) {
    std::cin >> i;
  }

  auto Merger = [](int i, int j) -> int { return i < j ? i : j; };
  SparseTable<int, static_cast<int (*)(int, int)>(Merger)> sparse_table(v);

  for (int i = 0; i < q; ++i) {
    int first;
    int last;
    std::cin >> first >> last;
    std::cout << sparse_table.Get(first, last) << "\n";
  }

  return 0;
}
