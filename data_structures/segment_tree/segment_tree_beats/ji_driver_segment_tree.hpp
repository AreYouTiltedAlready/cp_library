#include <vector>

namespace ds {
namespace segment_tree_beats {

template <typename T>
class JiDriverSegmentTree {
 public:
  struct Node {
    static constexpr T kInfinity = std::numeric_limits<T>::max() / 2;

    Node() : Node(0) {}
    explicit Node(T value)
        : max(value),
          min(value),
          sum(value),
          second_max(-kInfinity),
          second_min(kInfinity),
          lazy_add(0),
          count_max(1),
          count_min(1) {}

    friend Node operator+(const Node& lhs, const Node& rhs) {
      Node res{};
      res.sum = lhs.sum + rhs.sum;
      if (lhs.max == rhs.max) {
        res.max = lhs.max;
        res.count_max = lhs.count_max + rhs.count_max;
        res.second_max = std::max(lhs.second_max, rhs.second_max);
      } else if (lhs.max > rhs.max) {
        res.max = lhs.max;
        res.count_max = lhs.count_max;
        res.second_max = std::max(lhs.second_max, rhs.max);
      } else {
        res.max = rhs.max;
        res.count_max = rhs.count_max;
        res.second_max = std::max(lhs.max, rhs.second_max);
      }
      if (lhs.min == rhs.min) {
        res.min = lhs.min;
        res.count_min = lhs.count_min + rhs.count_min;
        res.second_min = std::min(lhs.second_min, rhs.second_min);
      } else if (lhs.min < rhs.min) {
        res.min = lhs.min;
        res.count_min = lhs.count_min;
        res.second_min = std::min(lhs.second_min, rhs.min);
      } else {
        res.min = rhs.min;
        res.count_min = rhs.count_min;
        res.second_min = std::min(lhs.min, rhs.second_min);
      }
      return res;
    }

    void ApplyAdd(T value, int left, int right) {
      max += value;
      min += value;
      lazy_add += value;
      sum += value * (right - left);
      if (second_max != -kInfinity) {
        second_max += value;
      }
      if (second_min != kInfinity) {
        second_min += value;
      }
    }

    void ApplyChmax(T value, int left, int right) {
      if (value <= min) {
        return;
      }
      sum += (value - min) * count_min;
      min = value;
      if (left + 1 == right) {
        max = value;
      } else {
        if (value >= max) {
          max = value;
        } else if (value > second_max) {
          second_max = value;
        }
      }
    }

    void ApplyChmin(T value, int left, int right) {
      if (value >= max) {
        return;
      }
      sum += (value - max) * count_max;
      max = value;
      if (left + 1 == right) {
        min = value;
      } else {
        if (value <= min) {
          min = value;
        } else if (value < second_min) {
          second_min = value;
        }
      }
    }

    T max;
    T min;
    T sum;
    T second_max;
    T second_min;
    T lazy_add;
    int count_max;
    int count_min;
  };

  explicit JiDriverSegmentTree(const std::vector<T>& values)
      : n_(static_cast<int>(values.size())), tree_(2 * n_) {
    Build(0, 0, n_, values);
  }

  void ApplyChmax(int left, int right, T value) {
    ApplyChmax(0, 0, n_, left, right, value);
  }

  void ApplyChmin(int left, int right, T value) {
    ApplyChmin(0, 0, n_, left, right, value);
  }

  void ApplyAdd(int left, int right, T value) {
    ApplyAdd(0, 0, n_, left, right, value);
  }

  Node Get(int left, int right) { return Get(0, 0, n_, left, right); }

 private:
  void Build(int x, int l, int r, const std::vector<T>& values) {
    if (l + 1 == r) {
      tree_[x] = Node(values[l]);
      return;
    }
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    Build(x + 1, l, mid, values);
    Build(z, mid, r, values);
    Pull(x, z);
  }

  Node Get(int x, int l, int r, int ql, int qr) {
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
    return Get(x + 1, l, mid, ql, qr) + Get(z, mid, r, ql, qr);
  }

  void ApplyChmax(int x, int l, int r, int ql, int qr, T value) {
    if (tree_[x].min >= value) {
      return;
    }
    if (ql <= l && r <= qr && tree_[x].second_min > value) {
      ApplyChmax(x, l, r, value);
      return;
    }
    Push(x, l, r);
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    if (ql < mid) {
      ApplyChmax(x + 1, l, mid, ql, qr, value);
    }
    if (qr > mid) {
      ApplyChmax(z, mid, r, ql, qr, value);
    }
    Pull(x, z);
  }

  void ApplyChmin(int x, int l, int r, int ql, int qr, T value) {
    if (tree_[x].max <= value) {
      return;
    }
    if (ql <= l && r <= qr && tree_[x].second_max < value) {
      ApplyChmin(x, l, r, value);
      return;
    }
    Push(x, l, r);
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    if (ql < mid) {
      ApplyChmin(x + 1, l, mid, ql, qr, value);
    }
    if (qr > mid) {
      ApplyChmin(z, mid, r, ql, qr, value);
    }
    Pull(x, z);
  }

  void ApplyAdd(int x, int l, int r, int ql, int qr, T value) {
    if (ql <= l && r <= qr) {
      ApplyAdd(x, l, r, value);
      return;
    }
    Push(x, l, r);
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    if (ql < mid) {
      ApplyAdd(x + 1, l, mid, ql, qr, value);
    }
    if (qr > mid) {
      ApplyAdd(z, mid, r, ql, qr, value);
    }
    Pull(x, z);
  }

  void ApplyChmax(int x, int l, int r, T value) {
    tree_[x].ApplyChmax(value, l, r);
  }

  void ApplyChmin(int x, int l, int r, T value) {
    tree_[x].ApplyChmin(value, l, r);
  }

  void ApplyAdd(int x, int l, int r, T value) {
    tree_[x].ApplyAdd(value, l, r);
  }

  void Push(int x, int l, int r) {
    int mid = (l + r) / 2;
    int z = x + (mid - l) * 2;
    if (T& add = tree_[x].lazy_add; add != 0) {
      ApplyAdd(x + 1, l, mid, add);
      ApplyAdd(z, mid, r, add);
      add = 0;
    }
    ApplyChmax(x + 1, l, mid, tree_[x].min);
    ApplyChmin(x + 1, l, mid, tree_[x].max);
    ApplyChmax(z, mid, r, tree_[x].min);
    ApplyChmin(z, mid, r, tree_[x].max);
  }

  inline void Pull(int x, int z) { tree_[x] = tree_[x + 1] + tree_[z]; }

  int n_;
  std::vector<Node> tree_;
};

}  // namespace segment_tree_beats
}  // namespace ds
