#include <limits>
#include <numeric>
#include <vector>

namespace ds {
namespace li_chao_tree {
// Tree for dealing with linear functions
// Supports line/segment insertion and min query at arbitrary point
// Note: it actually consumes an enormous amount of memory (basically, $O(q \log
// C)$, where C is the query range), but with segment insertion the bound
// increases up to $O(q {\log C}^2)$, so think twice about allocations
// `queries_count` in constructor has been added just as hint about approximate
// amount of memory
template <typename T>
class LiChaoTree {
 public:
  explicit LiChaoTree(T left, T right, int queries_count = 0)
      : lines_(), left_(left), right_(right) {
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
  class Line {
   public:
    explicit Line(T k = 0, T b = std::numeric_limits<T>::max()) : k(k), b(b) {}

    [[nodiscard]] T operator()(T x) const noexcept { return k * x + b; }

   private:
    T k;
    T b;
  };

  T GetMin(int x, T left, T right, T point) {
    T result = lines_[x](point);
    while (left + 1 != right) {
      const T middle = std::midpoint(left, right);
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
      const T middle = std::midpoint(left, right);
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
    const T middle = std::midpoint(left, right);
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
  int size_{};
};

}  // namespace li_chao_tree
}  // namespace ds
