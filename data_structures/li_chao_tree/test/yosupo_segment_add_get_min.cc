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

template <typename T>
class LiChaoTree {
 public:
  explicit LiChaoTree(T left, T right, T queries_count = 0)
      : lines_(),
        left_child_(),
        right_child_(),
        left_(left),
        right_(right),
        size_(0) {
    if (queries_count != 0) {
      const T lg = static_cast<T>(std::__lg(right - left + 1));
      lines_.reserve(queries_count * lg * lg);
      left_child_.reserve(queries_count * lg * lg);
      right_child_.reserve(queries_count * lg * lg);
    }
    MakeNode();
    MakeNode();
  }

  T GetMin(T point) { return GetMin(1, left_, right_, point); }

  void InsertSegment(T left, T right, T k, T b) {
    InsertSegment(1, left_, right_, left, right, Line(k, b));
  }

  void InsertLine(T k, T b) { InsertLine(1, left_, right_, Line(k, b)); }

 private:
  struct Line {
    Line(T k = 0, T b = std::numeric_limits<T>::max()) : k(k), b(b) {}

    [[nodiscard]] T operator()(T x) const noexcept { return k * x + b; }

    T k;
    T b;
  };

  T GetMin(int x, T left, T right, T point) {
    T result = lines_[x](point);
    while (left + 1 != right) {
      const T middle = left + (right - left) / 2;
      if (point < middle) {
        ExtendLeft(x);
        x = left_child_[x];
        right = middle;
      } else {
        ExtendRight(x);
        x = right_child_[x];
        left = middle;
      }
      result = std::min(result, lines_[x](point));
    }
    return result;
  }

  void InsertLine(int x, T left, T right, Line line) {
    while (left + 1 != right) {
      const T middle = left + (right - left) / 2;
      const bool new_less_on_middle = line(middle) < lines_[x](middle);
      const bool new_less_on_left = line(left) < lines_[x](left);
      if (new_less_on_middle) {
        std::swap(lines_[x], line);
      }

      if (new_less_on_left == new_less_on_middle) {
        ExtendRight(x);
        x = right_child_[x];
        left = middle;
      } else {
        ExtendLeft(x);
        x = left_child_[x];
        right = middle;
      }
    }

    if (lines_[x](left) > line(left)) {
      std::swap(lines_[x], line);
    }
  }

  void InsertSegment(int x, T left, T right, T segment_left, T segment_right,
                     Line line) {
    if (segment_left <= left && right <= segment_right) {
      InsertLine(x, left, right, line);
      return;
    }
    const T middle = left + (right - left) / 2;
    if (segment_left < middle) {
      ExtendLeft(x);
      InsertSegment(left_child_[x], left, middle, segment_left, segment_right,
                    line);
    }
    if (segment_right > middle) {
      ExtendRight(x);
      InsertSegment(right_child_[x], middle, right, segment_left, segment_right,
                    line);
    }
  }

  void ExtendLeft(int x) {
    if (left_child_[x] != 0) {
      return;
    }
    left_child_[x] = MakeNode();
  }

  void ExtendRight(int x) {
    if (right_child_[x] != 0) {
      return;
    }
    right_child_[x] = MakeNode();
  }

  int MakeNode() {
    lines_.emplace_back();
    left_child_.push_back(0);
    right_child_.push_back(0);
    return size_++;
  }

  std::vector<Line> lines_;
  std::vector<int> left_child_;
  std::vector<int> right_child_;

  T left_;
  T right_;
  int size_;
};

// https://judge.yosupo.jp/problem/segment_add_get_min
// https://judge.yosupo.jp/submission/180782
void RunCase([[maybe_unused]] int testcase) {
  int n;
  int q;
  std::cin >> n >> q;

  constexpr int64_t kRange = 1e9 + 5;
  LiChaoTree<int64_t> li_chao_tree(-kRange, kRange, n + q);
  for (int i = 0; i < n; ++i) {
    int l;
    int r;
    int64_t a;
    int64_t b;
    std::cin >> l >> r >> a >> b;
    li_chao_tree.InsertSegment(l, r, a, b);
  }

  while (q--) {
    int t;
    std::cin >> t;
    if (t == 0) {
      int l;
      int r;
      int64_t k;
      int64_t b;
      std::cin >> l >> r >> k >> b;
      li_chao_tree.InsertSegment(l, r, k, b);
    } else {
      int64_t p;
      std::cin >> p;
      if (int64_t result = li_chao_tree.GetMin(p);
          result == std::numeric_limits<int64_t>::max()) {
        std::cout << "INFINITY\n";
      } else {
        std::cout << result << "\n";
      }
    }
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
