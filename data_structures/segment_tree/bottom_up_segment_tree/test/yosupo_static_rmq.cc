#include <bits/stdc++.h>

namespace {

template <class Fun>
class y_combinator_result {
  Fun fun_;

 public:
  template <class T>
  explicit y_combinator_result(  // NOLINT(*-forwarding-reference-overload)
      T&& fun)
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

// Simple bottom-up segment tree
// Useful for easy invariants only
template <typename S, auto Op, auto E>
class SegmentTree {
  static_assert(std::is_convertible_v<decltype(Op), S (*)(S, S)>,
                "Op must work as S(S, S)");
  static_assert(std::is_convertible_v<decltype(E), S (*)()>,
                "E must work as S()");

 public:
  explicit SegmentTree(int n) : n_(n), tree_(n * 2, E()) {}

  template <typename U, typename std::enable_if_t<std::is_assignable_v<S&, U>,
                                                  void>* = nullptr>
  explicit SegmentTree(const std::vector<U>& values)
      : SegmentTree(static_cast<int>(values.size())) {
    std::copy(values.cbegin(), values.cend(), tree_.begin() + n_);
    for (int i = n_ - 1; i > 0; --i) {
      Pull(i);
    }
  }

  void Apply(int pos, const S& value) {
    tree_[pos += n_] = value;
    pos /= 2;
    while (pos > 0) {
      Pull(pos);
      pos /= 2;
    }
  }

  S Get(int left, int right) const {
    left += n_;
    right += n_;
    S res = E();
    while (left < right) {
      if (left % 2 == 1) {
        res = Op(res, tree_[left++]);
      }
      if (right % 2 == 1) {
        res = Op(res, tree_[--right]);
      }
      left /= 2;
      right /= 2;
    }
    return res;
  }

 private:
  inline void Pull(int k) { tree_[k] = Op(tree_[k * 2], tree_[k * 2 + 1]); }

  int n_;
  std::vector<S> tree_;
};

inline int Min(int lhs, int rhs) { return lhs < rhs ? lhs : rhs; }
inline int E() { return std::numeric_limits<int>::max(); }

// https://judge.yosupo.jp/problem/staticrmq
// https://judge.yosupo.jp/submission/179660
void RunCase([[maybe_unused]] int testcase) {
  int n;
  int q;
  std::cin >> n >> q;

  std::vector<int> v(n);
  for (int& i : v) {
    std::cin >> i;
  }

  SegmentTree<int, Min, E> segment_tree(v);
  while (q--) {
    int left;
    int right;
    std::cin >> left >> right;
    std::cout << segment_tree.Get(left, right) << "\n";
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
