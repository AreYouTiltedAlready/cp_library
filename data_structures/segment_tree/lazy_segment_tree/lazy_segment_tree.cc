#include <functional>
#include <vector>

// Fully generic implementation of lazy segment tree
// Unlike bottom-up implementation, this one does not need neutral element
// * $S$ - some semigroup
// * $T$ - lazy tag which we want to apply. Tag() MUST yield the tag which maps
// any item from S to itself (e.g. 0 for addition, 1 for multiplication, ...)
// * $Op$ - $S x S -> S$ - associative binary operation on
// given semigroup
// * $Mapping$ - $ S x T -> S (last two arguments may be int, int - check range
// affine range sum for better understanding)$ - function which defines
// behaviour of tag (addition on segment, affine etc.). The last two ints stand
// for range bounds (which helps to save memory not storing the length
// internally in semigroup S)
// * Composition - $T x T -> T$ defines the composition rules of two tags in
// order they are applied
// FindFirst(first, last, P) and FindLast(first, last, P)
// find first and last elements (correspondingly) in range on which P yields
// true.
// $P$ MUST be monotonic (it means, if $P$ is true for some $(L, R)$, then P
// must be true for $(L', R')$ : $L' <= L < R <= R'$
// They can also be used for
// searching smth like lower_bound For example, we store in $S$ the sum on
// corresponding range, we have range $(L, R)$ and we want to find first (or
// last) index $i$, s.t. Sum[L, i] = k Then we can modify P in the following
// way: function P(S):
//  if Sum(S) < k:
//    k -= Sum(s)
//    return false // Go to another direction
//  return true // We will go deeper and find the sought position
// Everything works in $O(\log n)$
template <typename S, typename T, auto Op, auto Mapping, auto Composition>
class LazySegmentTree {
  static_assert(
      std::is_convertible_v<decltype(Op), S (*)(S, S)> ||
          std::is_convertible_v<decltype(Op), S (*)(const S&, const S&)>,
      "Op must work as S(S, S)");
  static_assert(
      std::is_convertible_v<decltype(Mapping), S (*)(S, T, int, int)> ||
          std::is_convertible_v<decltype(Mapping),
                                S (*)(const S&, const T&, int, int)>,
      "Mapping must work as S(S, S, int(?), int(?)");
  static_assert(std::is_convertible_v<decltype(Composition), T (*)(T, T)> ||
                    std::is_convertible_v<decltype(Composition),
                                          T (*)(const T&, const T&)>,
                "Composition must work as T(T, T)");

 public:
  explicit LazySegmentTree(int n) : n_(n), tree_(n * 2), tag_(n * 2) {}

  template <typename U, typename std::enable_if_t<std::is_assignable_v<S&, U>,
                                                  void>* = nullptr>
  explicit LazySegmentTree(const std::vector<U>& values)
      : LazySegmentTree(static_cast<int>(values.size())) {
    Build(0, 0, n_, values);
  }

  template <typename U, typename std::enable_if_t<std::is_assignable_v<S&, U>,
                                                  void>* = nullptr>
  void Apply(int pos, const U& value) {
    Apply(0, 0, n_, pos, value);
  }

  void Apply(int left, int right, const T& tag) {
    Apply(0, 0, n_, left, right, tag);
  }

  [[nodiscard]] S Get(int pos) { return Get(0, 0, n_, pos); }

  [[nodiscard]] S Get(int left, int right) {
    return Get(0, 0, n_, left, right);
  }

  [[nodiscard]] S GetAll() const { return tree_[0]; }

  using Predicate = std::function<bool(const S&)>;

  int FindFirst(int left, int right, const Predicate& pred) {
    return FindFirst(0, 0, n_, left, right, pred);
  }

  int FindLast(int left, int right, const Predicate& pred) {
    return FindLast(0, 0, n_, left, right, pred);
  }

 private:
  template <typename U, typename std::enable_if_t<std::is_assignable_v<S&, U>,
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

  template <typename U, typename std::enable_if_t<std::is_assignable_v<S&, U>,
                                                  void>* = nullptr>
  void Apply(int x, int l, int r, int pos, const U& value) {
    if (l + 1 == r) {
      tree_[x] = value;
      return;
    }
    Push(x, l, r);
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    if (pos < mid) {
      Apply(x + 1, l, mid, pos, value);
    } else {
      Apply(z, mid, r, pos, value);
    }
    Pull(x, z);
  }

  void Apply(int x, int l, int r, int ql, int qr, const T& tag) {
    if (ql <= l && r <= qr) {
      Apply(x, l, r, tag);
      return;
    }
    Push(x, l, r);
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    if (ql < mid) {
      Apply(x + 1, l, mid, ql, qr, tag);
    }
    if (qr > mid) {
      Apply(z, mid, r, ql, qr, tag);
    }
    Pull(x, z);
  }

  [[nodiscard]] S Get(int x, int l, int r, int pos) {
    while (l + 1 != r) {
      Push(x, l, r);
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

  [[nodiscard]] S Get(int x, int l, int r, int ql, int qr) {
    if (ql <= l && r <= qr) {
      return tree_[x];
    }
    Push(x, l, r);
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
    Push(x, l, r);
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
    Push(x, l, r);
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
      Push(x, l, r);
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
      Push(x, l, r);
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

  void Apply(int x, int l, int r, const T& tag) {
    tree_[x] = Mapping(tree_[x], tag, l, r);
    tag_[x] = Composition(tag_[x], tag);
  }

  void Push(int x, int l, int r) {
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    Apply(x + 1, l, mid, tag_[x]);
    Apply(z, mid, r, tag_[x]);
    tag_[x] = {};
  }

  inline void Pull(int x, int z) { tree_[x] = Op(tree_[x + 1], tree_[z]); }

  int n_;
  std::vector<T> tag_;
  std::vector<S> tree_;
};
