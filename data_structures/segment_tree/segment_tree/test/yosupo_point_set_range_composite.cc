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

template <typename T, auto Op>
class SegmentTree {
  static_assert(
      std::is_convertible_v<decltype(Op), T (*)(T, T)> ||
          std::is_convertible_v<decltype(Op), T (*)(const T&, const T&)>,
      "Op must work as T(T, T)");

 public:
  explicit SegmentTree(int n) : n_(n), tree_(n * 2) {}

  template <typename U, typename std::enable_if_t<std::is_assignable_v<T&, U>,
                                                  void>* = nullptr>
  explicit SegmentTree(const std::vector<U>& values)
      : SegmentTree(static_cast<int>(values.size())) {
    Build(0, 0, n_, values);
  }

  template <typename U, typename std::enable_if_t<std::is_assignable_v<T&, U>,
                                                  void>* = nullptr>
  void Apply(int pos, const U& value) {
    Apply(0, 0, n_, pos, value);
  }

  [[nodiscard]] T Get(int pos) const { return Get(0, 0, n_, pos); }

  [[nodiscard]] T Get(int left, int right) const {
    return Get(0, 0, n_, left, right);
  }

  [[nodiscard]] T GetAll() const { return tree_[0]; }

  using Predicate = std::function<bool(const T&)>;

  int FindFirst(int left, int right, const Predicate& pred) {
    return FindFirst(0, 0, n_, left, right, pred);
  }

  int FindLast(int left, int right, const Predicate& pred) {
    return FindLast(0, 0, n_, left, right, pred);
  }

 private:
  template <typename U, typename std::enable_if_t<std::is_assignable_v<T&, U>,
                                                  void>* = nullptr>
  void Build(int x, int l, int r, const std::vector<U>& values) {
    if (l + 1 == r) {
      tree_[x] = values[l];
      return;
    }
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    Build(x + 1, l, mid, values);
    Build(z, mid, r, values);
    Pull(x, z);
  }

  template <typename U, typename std::enable_if_t<std::is_assignable_v<T&, U>,
                                                  void>* = nullptr>
  void Apply(int x, int l, int r, int pos, const U& value) {
    if (l + 1 == r) {
      tree_[x] = value;
      return;
    }
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    if (pos < mid) {
      Apply(x + 1, l, mid, pos, value);
    } else {
      Apply(z, mid, r, pos, value);
    }
    Pull(x, z);
  }

  [[nodiscard]] T Get(int x, int l, int r, int pos) const {
    while (l + 1 != r) {
      int mid = (l + r) / 2;
      int z = x + (mid - l) * 2;
      if (pos < mid) {
        r = mid;
        x = x + 1;
      } else {
        l = mid;
        x = z;
      }
    }
    return tree_[x];
  }

  [[nodiscard]] T Get(int x, int l, int r, int ql, int qr) const {
    if (ql <= l && r <= qr) { return tree_[x]; }
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    if (qr <= mid) { return Get(x + 1, l, mid, ql, qr); }
    if (ql >= mid) { return Get(z, mid, r, ql, qr); }
    return Op(Get(x + 1, l, mid, ql, qr), Get(z, mid, r, ql, qr));
  }

  int FindFirst(int x, int l, int r, int ql, int qr, const Predicate& pred) {
    if (ql <= l && r <= qr) {
      if (!pred(tree_[x])) { return -1; }
      return FindFirstKnowingly(x, l, r, pred);
    }
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    int result = -1;
    if (ql < mid) { result = FindFirst(x + 1, l, mid, ql, qr, pred); }
    if (result == -1 && qr > mid) {
      result = FindFirst(z, mid, r, ql, qr, pred);
    }
    return result;
  }

  int FindLast(int x, int l, int r, int ql, int qr, const Predicate& pred) {
    if (ql <= l && r <= qr) {
      if (!pred(tree_[x])) { return -1; }
      return FindLastKnowingly(x, l, r, pred);
    }
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    int result = -1;
    if (qr > mid) { result = FindLast(z, mid, r, ql, qr, pred); }
    if (result == -1 && ql < mid) {
      result = FindLast(x + 1, l, mid, ql, qr, pred);
    }
    return result;
  }

  int FindFirstKnowingly(int x, int l, int r, const Predicate& pred) {
    while (l + 1 != r) {
      int mid = (l + r) / 2;
      int z = x + (mid - l) * 2;
      if (pred(tree_[x + 1])) {
        r = mid;
        x = x + 1;
      } else {
        l = mid;
        x = z;
      }
    }
    return l;
  }

  int FindLastKnowingly(int x, int l, int r, const Predicate& pred) {
    while (l + 1 != r) {
      int mid = (l + r) / 2;
      int z = x + (mid - l) * 2;
      if (pred(tree_[z])) {
        l = mid;
        x = z;
      } else {
        r = mid;
        x = x + 1;
      }
    }
    return l;
  }

  inline void Pull(int x, int z) { tree_[x] = Op(tree_[x + 1], tree_[z]); }

  int n_;
  std::vector<T> tree_;
};

class Barrett {
 public:
  constexpr explicit Barrett(uint32_t mod)
      : mod_inverse(static_cast<uint64_t>(-1) / mod + 1), mod(mod) {}

  [[nodiscard]] uint32_t Product(uint32_t lhs, uint32_t rhs) const {
    return (*this)(static_cast<uint64_t>(lhs) * rhs);
  }

  [[nodiscard]] uint32_t operator()(uint64_t n) const {
    auto x = static_cast<uint64_t>(
        (static_cast<__uint128_t>(n) * mod_inverse) >> 64);
    uint64_t m = x * mod;
    return n - m + (n < m ? mod : 0);
  }

 private:
  uint64_t mod_inverse;
  uint32_t mod;
};

template <uint32_t kMod>
class StaticModint {
 public:
  static_assert(kMod < (1U << 30));

  StaticModint() : value_(0) {}

  template <typename T,
            typename std::enable_if_t<
                std::is_integral_v<T> && std::is_signed_v<T>, void>* = nullptr>
  StaticModint(T value)
      : value_(  // NOLINT(*-explicit-constructor)
            value % kMod) {
    if (value < 0) { value += kMod; }
  }

  template <typename T, typename std::enable_if_t<std::is_integral_v<T> &&
                                                      std::is_unsigned_v<T>,
                                                  void>* = nullptr>
  StaticModint(T value)
      : value_(  // NOLINT(*-explicit-constructor)
            value % kMod) {}

  using Mint = StaticModint<kMod>;

  static Mint Raw(unsigned value) {
    Mint result;
    result.value_ = value;
    return result;
  }

  [[nodiscard]] Mint operator+() const noexcept { return Mint(*this); }

  [[nodiscard]] Mint operator-() const noexcept { return Mint(kMod - value_); }

  [[nodiscard]] bool IsZero() const noexcept { return value_ == 0; }

  [[nodiscard]] explicit operator uint32_t() const { return value_; }

  [[nodiscard]] Mint Inverse() const noexcept {
    Mint m(*this);
    Mint result = Raw(1U);
#pragma GCC unroll(30)
    for (unsigned i = 0; i < 30; ++i) {
      if (((kMod - 2) >> i) % 2 == 1) { result *= m; }
      m *= m;
    }
    return result;
  }

  Mint& operator+=(const Mint& other) {
    if (value_ += other.value_; value_ >= kMod) { value_ -= kMod; }
    return *this;
  }

  Mint& operator-=(const Mint& other) {
    if (value_ += kMod - other.value_; value_ >= kMod) { value_ -= kMod; }
    return *this;
  }

  Mint& operator*=(const Mint& other) {
    value_ = barrett.Product(value_, other.value_);
    return *this;
  }

  Mint& operator/=(const Mint& other) {
    value_ = barrett.Product(value_, other.Inverse().value_);
    return *this;
  }

  friend Mint operator+(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) += rhs;
  }

  friend Mint operator-(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) -= rhs;
  }

  friend Mint operator*(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) *= rhs;
  }

  friend Mint operator/(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) /= rhs;
  }

  friend Mint Power(Mint m, uint64_t n) {
    Mint result = Raw(1U);
    while (n > 0) {
      if (n % 2 == 1) { result *= m; }
      m *= m;
      n /= 2;
    }
    return result;
  }

  friend bool operator==(const Mint& lhs, const Mint& rhs) {
    return lhs.value_ == rhs.value_;
  }

  friend std::istream& operator>>(std::istream& is, Mint& mint) {
    uint32_t value;
    is >> value;
    mint = Mint(value);
    return is;
  }

  friend std::ostream& operator<<(std::ostream& os, const Mint& mint) {
    return os << mint.value_;
  }

 private:
  static constexpr Barrett barrett = Barrett(kMod);
  uint32_t value_;
};

using Mint = StaticModint<998244353>;

struct S {
  S() : a(), b() {}
  explicit S(Mint a, Mint b) : a(a), b(b) {}

  [[nodiscard]] Mint operator()(Mint x) const { return a * x + b; }

  friend std::istream& operator>>(std::istream& is, S& s) {
    return is >> s.a >> s.b;
  }

  Mint a;
  Mint b;
};

S Op(const S& lhs, const S& rhs) {
  return S(lhs.a * rhs.a, lhs.b * rhs.a + rhs.b);
}

// https://judge.yosupo.jp/problem/point_set_range_composite
// https://judge.yosupo.jp/submission/179417
void RunCase([[maybe_unused]] int testcase) {
  int n;
  int q;
  std::cin >> n >> q;

  std::vector<S> data(n);
  for (S& s : data) { std::cin >> s; }

  SegmentTree<S, Op> segment_tree(data);
  while (q--) {
    int t;
    std::cin >> t;
    if (t == 0) {
      int p;
      Mint x;
      Mint y;
      std::cin >> p >> x >> y;
      segment_tree.Apply(p, S(x, y));
    } else {
      int left;
      int right;
      Mint x;
      std::cin >> left >> right >> x;
      std::cout << segment_tree.Get(left, right)(x) << "\n";
    }
  }
}

void Main() {
  int testcases = 1;
  // std::cin >> testcases;
  for (int tt = 1; tt <= testcases; ++tt) { RunCase(tt); }
}

}  // namespace

int main() {
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);
  Main();
  return 0;
}
