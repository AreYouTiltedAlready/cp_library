#include <functional>
#include <vector>

template <typename T, T (*Op)(T, T)>
class SegmentTree {
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
    ;
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
    if (ql <= l && r <= qr) {
      return tree_[x];
    }
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    if (qr <= mid) {
      return Get(x + 1, l, mid, ql, qr);
    }
    if (ql >= mid) {
      return Get(z, mid, r, ql, qr);
    }
    return Op(Get(x + 1, l, mid, ql, qr), Get(z, mid, r, ql, qr));
  }

  int FindFirst(int x, int l, int r, int ql, int qr, const Predicate& pred) {
    if (ql <= l && r <= qr) {
      if (!pred(tree_[x])) {
        return -1;
      }
      return FindFirstKnowingly(x, l, r, pred);
    }
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    int result = -1;
    if (ql < mid) {
      result = FindFirst(x + 1, l, mid, ql, qr, pred);
    }
    if (result == -1 && qr > mid) {
      result = FindFirst(z, mid, r, ql, qr, pred);
    }
    return result;
  }

  int FindLast(int x, int l, int r, int ql, int qr, const Predicate& pred) {
    if (ql <= l && r <= qr) {
      if (!pred(tree_[x])) {
        return -1;
      }
      return FindLastKnowingly(x, l, r, pred);
    }
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    int result = -1;
    if (qr > mid) {
      result = FindLast(z, mid, r, ql, qr, pred);
    }
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
