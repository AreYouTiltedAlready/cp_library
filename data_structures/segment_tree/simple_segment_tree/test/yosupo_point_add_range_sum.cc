// Problem: https://judge.yosupo.jp/problem/point_add_range_sum
// Submission: https://judge.yosupo.jp/submission/182093

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
namespace segment_tree {

namespace internal {

template <typename T>
constexpr bool has_binary_plus = requires(const T& lhs, const T& rhs) {
  { lhs + rhs } -> std::same_as<T>;
};

template <typename T>
concept monoid = std::is_default_constructible_v<T> && has_binary_plus<T>;

}  // namespace internal

template <internal::monoid S>
class SegmentTree {
 public:
  explicit SegmentTree(int n) : tree_(n * 2), n_(n) {}

  template <typename U>
  requires std::is_assignable_v<S&, U>
  explicit SegmentTree(int n, const U& value) : SegmentTree(n) {
    std::ranges::fill(tree_, value);
    for (int i = n_ - 1; i > 0; --i) {
      Pull(i);
    }
  }

  template <typename U>
  requires std::is_assignable_v<S&, U>
  explicit SegmentTree(const std::vector<U>& values)
      : SegmentTree(static_cast<int>(values.size())) {
    std::ranges::copy(values, tree_.begin() + n_);
    for (int i = n_ - 1; i > 0; --i) {
      Pull(i);
    }
  }

  template <typename U>
  requires std::is_assignable_v<S&, U>
  void Set(int pos, const U& value) {
    pos += n_;
    tree_[pos] = value;
    DoPull(pos);
  }

  void Apply(int pos, const S& value) {
    pos += n_;
    tree_[pos] = tree_[pos] + value;
    DoPull(pos);
  }

  [[nodiscard]] S Get() const { return tree_[1]; }

  [[nodiscard]] S Get(int pos) const {
    pos += n_;
    return tree_[pos];
  }

  [[nodiscard]] S Get(int first, int last) const {
    first += n_;
    last += n_;
    S res{};
    while (first < last) {
      if ((first & 1) == 1) {
        res = res + tree_[first++];
      }
      if ((last & 1) == 1) {
        res = res + tree_[--last];
      }
      first >>= 1;
      last >>= 1;
    }

    return res;
  }

 private:
  void DoPull(int v) {
    v /= 2;
    while (v > 0) {
      Pull(v);
      v /= 2;
    }
  }

  inline void Pull(int k) { tree_[k] = tree_[k << 1] + tree_[k << 1 | 1]; }

  std::vector<S> tree_;
  int n_;
};

}  // namespace segment_tree

}  // namespace ds


void RunCase([[maybe_unused]] int testcase) {
  int n{};
  int q{};
  std::cin >> n >> q;

  std::vector<int> v(n);
  for (int& it : v) {
    std::cin >> it;
  }

  ds::segment_tree::SegmentTree<int64_t> fenwick(v);
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
