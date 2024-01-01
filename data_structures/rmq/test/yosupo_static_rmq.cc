// https://judge.yosupo.jp/problem/staticrmq
// https://judge.yosupo.jp/submission/180871

#include <algorithm>
#include <bit>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

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

// Range M(ax)|(in)imum Query problem solver
// Time: $O(n)$ / $O(1)$
// GetIndex(first, last) always yields the FIRST occurrence of min/max in range
namespace rmq {

namespace internal {

enum class RMQMode {
  kMax,
  kMin,
};

template <typename T, RMQMode mode>
class RMQSolver {
 public:
  static constexpr uint32_t kBlockLength = 32;
  RMQSolver()
      : values_(0), small_blocks_(), blocks_table_(), n_(0), blocks_count_(0) {}

  template <typename U,
            std::enable_if_t<std::is_same_v<std::decay_t<U>, std::vector<T>>,
                             void>* = nullptr>
  explicit RMQSolver(U&& values)
      : values_(std::forward<U>(values)),
        small_blocks_(static_cast<int>(values_.size())),
        n_(static_cast<int>(values_.size())),
        blocks_count_((static_cast<int>(values_.size()) + kBlockLength - 1) /
                      kBlockLength) {
    // building last 32 elements of monotonic stack for each index
    {
      std::vector<int> stack;
      stack.reserve(n_);
      stack.push_back(0);
      small_blocks_[0] = 1;
      for (int i = 1; i < n_; ++i) {
        small_blocks_[i] = small_blocks_[i - 1] << 1;
        while (!stack.empty() && values_[stack.back()] >= values_[i]) {
          if (int d = i - stack.back(); d < kBlockLength) {
            small_blocks_[i] ^= 1U << d;
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
        blocks_table_[0][i] = static_cast<int>(i * kBlockLength);
      }
      for (int i = 0; i < n_; ++i) {
        blocks_table_[0][i / kBlockLength] =
            Merger(blocks_table_[0][i / kBlockLength], i);
      }
      for (int i = 1; i < blocks_log; ++i) {
        for (int j = 0; j < blocks_count_ - (1 << i) + 1; ++j) {
          blocks_table_[i][j] =
              Merger(blocks_table_[i - 1][j],
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
      result = Merger(result, GetSmallBlock(first + kBlockLength - 1 - first_r,
                                            kBlockLength - first_r));
      first_q += 1;
    }
    if (last_r != 0) {
      result = Merger(result, GetSmallBlock(last - 1, last_r));
    }
    if (first_q < last_q) {
      result = Merger(result, GetOnBlocks(first_q, last_q));
    }

    return result;
  }

  [[nodiscard]] T GetValue(int first, int last) const {
    return values_[GetIndex(first, last)];
  }

 private:
  [[nodiscard]] inline int GetSmallBlock(int right, int length) const {
    return right + 1 -
           std::bit_width(small_blocks_[right] & ((1U << length) - 1));
  }

  [[nodiscard]] inline int GetOnBlocks(int first, int last) const {
    int level = std::bit_width(static_cast<uint32_t>(last - first)) - 1;
    return Merger(blocks_table_[level][first],
                  blocks_table_[level][last - (1 << level)]);
  }

  [[nodiscard]] inline int Merger(int lhs, int rhs) const {
    if constexpr (mode == RMQMode::kMin) {
      return values_[lhs] < values_[rhs] ? lhs : rhs;
    }
    return values_[lhs] > values_[rhs] ? lhs : rhs;
  }

  std::vector<T> values_;
  std::vector<uint32_t> small_blocks_;
  std::vector<std::vector<int>> blocks_table_;

  int n_;
  int blocks_count_;
};

}  // namespace internal

template <typename T>
using RMQMaxSolver = internal::RMQSolver<T, internal::RMQMode::kMax>;

template <typename T>
using RMQMinSolver = internal::RMQSolver<T, internal::RMQMode::kMin>;

}  // namespace rmq

void RunCase([[maybe_unused]] int testcase) {
  int n;
  int q;
  std::cin >> n >> q;

  std::vector<int> values(n);
  for (int& it : values) {
    std::cin >> it;
  }

  rmq::RMQMinSolver<int> rmq_solver(std::move(values));
  while (q--) {
    int first;
    int last;
    std::cin >> first >> last;
    std::cout << rmq_solver.GetValue(first, last) << "\n";
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
