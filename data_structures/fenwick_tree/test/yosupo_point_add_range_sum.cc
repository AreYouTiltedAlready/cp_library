// Problem: https://judge.yosupo.jp/problem/point_add_range_sum
// Submission: https://judge.yosupo.jp/submission/182092

#include <bits/stdc++.h>

namespace {

template <class Fun>
class y_combinator_result {
  Fun fun_;

 public:
  template <class T>
  explicit y_combinator_result(T&& fun)  // NOLINT(*forwarding-reference*)
      : fun_(std::forward<T>(fun)) {}

  template <class... Args>
  decltype(auto) operator()(Args&&... args) {
    return fun_(std::ref(*this), std::forward<Args>(args)...);
  }
};

template <class Fun>
decltype(auto) y_combinator(Fun&& fun) {
  return y_combinator_result<std::decay_t<Fun>>(std::forward<Fun>(fun));
}

namespace ds {
namespace fenwick_tree {

template <typename S>
class FenwickTree {
 public:
  explicit FenwickTree(int n) : tree_(n + 1), n_(n) {}

  void Apply(int pos, const S& value) {
    pos += 1;
    while (pos <= n_) {
      tree_[pos] += value;
      pos += pos & -pos;
    }
  }

  [[nodiscard]] S Get(int last) const {
    S res{};
    while (last > 0) {
      res += tree_[last];
      last -= last & -last;
    }
    return res;
  }

  [[nodiscard]] S Get(int first, int last) const {
    return Get(last) - Get(first);
  }

  [[nodiscard]] int LowerBound(S value) const {
    int res = 0;
    for (int i = std::bit_width(static_cast<uint32_t>(n_)); i >= 0; --i) {
      if (int next = res + (1 << i); next <= n_ && tree_[next] < value) {
        res = next;
        value -= tree_[res];
      }
    }
    return res;
  }

 private:
  std::vector<S> tree_;
  int n_;
};

}  // namespace fenwick_tree
}  // namespace ds

void RunCase([[maybe_unused]] int testcase) {
  int n{};
  int q{};
  std::cin >> n >> q;

  std::vector<int> v(n);
  for (int& it : v) {
    std::cin >> it;
  }

  ds::fenwick_tree::FenwickTree<int64_t> fenwick(n);
  for (int i = 0; i < n; ++i) {
    fenwick.Apply(i, v[i]);
  }

  for (int i = 0; i < q; ++i) {
    int t{};
    int a{};
    int b{};
    std::cin >> t >> a >> b;
    if (t == 0) {
      fenwick.Apply(a, b);
    } else {
      std::cout << fenwick.Get(a, b) << "\n";
    }
  }
}

void Main() {
  int testcases = 1;
  //std::cin >> testcases;
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
