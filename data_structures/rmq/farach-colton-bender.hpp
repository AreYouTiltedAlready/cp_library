#include <cstdint>
#include <cstdlib>
#include <vector>

namespace ds {
namespace rmq {
namespace internal {

template <typename T, typename Comp>
class RMQSolver {
 public:
  static constexpr int kBlockLength = 64;

  template <typename U,
            std::enable_if_t<std::is_same_v<std::decay_t<U>, std::vector<T>>,
                             void>* = nullptr>
  explicit RMQSolver(U&& values)
      : values_(std::forward<U>(values)),
        small_blocks_(values_.size()),
        comp_(),
        n_(static_cast<int>(values_.size())),
        blocks_count_((values_.size() + kBlockLength - 1) / kBlockLength) {
    // building last 64 elements of monotonic stack for each index
    {
      std::vector<int> stack;
      stack.reserve(n_);
      stack.push_back(0);
      small_blocks_[0] = 1;
      for (int i = 1; i < n_; ++i) {
        small_blocks_[i] = small_blocks_[i - 1] << 1;
        while (!stack.empty() && values_[stack.back()] >= values_[i]) {
          if (int d = i - stack.back(); d < kBlockLength) {
            small_blocks_[i] ^= static_cast<uint64_t>(1) << d;
          }
          stack.pop_back();
        }
        stack.push_back(i);
        small_blocks_[i] += 1;
      }
    }

    // Building the sparse table for blocks of size n / kBlockLength
    {
      const int blocks_log =
          std::bit_width(static_cast<uint32_t>(blocks_count_));
      blocks_table_.resize(blocks_log);
      for (int i = 0; i < blocks_log; ++i) {
        blocks_table_[i].resize(blocks_count_ - (1 << i) + 1);
      }
      for (int i = 0; i < blocks_count_; ++i) {
        blocks_table_[0][i] = i * kBlockLength;
      }
      for (int i = 0; i < n_; ++i) {
        blocks_table_[0][i / kBlockLength] =
            Merge(blocks_table_[0][i / kBlockLength], i);
      }
      for (int i = 1; i < blocks_log; ++i) {
        for (int j = 0; j < blocks_count_ - (1 << i) + 1; ++j) {
          blocks_table_[i][j] = Merge(blocks_table_[i - 1][j],
                                      blocks_table_[i - 1][j + (1 << (i - 1))]);
        }
      }
    }
  }

  [[nodiscard]] int GetIndex(int first, int last) const {
    auto [first_q, first_r] = std::div(first, kBlockLength);
    auto [last_q, last_r] = std::div(last, kBlockLength);

    if (first_q == last_q) {
      return GetSmallBlock(last - 1, last - first);
    }

    int result = first;
    if (first_r != 0) {
      result = Merge(result, GetSmallBlock(first + kBlockLength - 1 - first_r,
                                            kBlockLength - first_r));
      first_q += 1;
    }
    if (last_r != 0) {
      result = Merge(result, GetSmallBlock(last - 1, last_r));
    }
    if (first_q < last_q) {
      result = Merge(result, GetOnBlocks(first_q, last_q));
    }

    return result;
  }

  [[nodiscard]] T GetValue(int first, int last) const {
    return values_[GetIndex(first, last)];
  }

 private:
  [[nodiscard]] inline int GetSmallBlock(int right, int length) const {
    return right + 1 -
           std::bit_width(small_blocks_[right] &
                          ((static_cast<uint64_t>(1) << length) - 1));
  }

  [[nodiscard]] inline int GetOnBlocks(int first, int last) const {
    int level = std::bit_width(static_cast<uint32_t>(last - first)) - 1;
    return Merge(blocks_table_[level][first],
                 blocks_table_[level][last - (1 << level)]);
  }

  [[nodiscard]] inline int Merge(int lhs, int rhs) const {
    return comp_(values_[lhs], values_[rhs]) ? lhs : rhs;
  }

  std::vector<T> values_;
  std::vector<uint64_t> small_blocks_;
  std::vector<std::vector<int>> blocks_table_;

  const Comp comp_;

  int n_;
  int blocks_count_;
};

}  // namespace internal

template <typename T>
using RMQMinSolver = internal::RMQSolver<T, std::less<>>;

template <typename T>
using RMQMaxSolver = internal::RMQSolver<T, std::greater<>>;

}  // namespace rmq
}  // namespace ds
