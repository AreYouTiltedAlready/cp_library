#include <algorithm>
#include <numeric>
#include <string_view>
#include <vector>

class SuffixArray {
 public:
  explicit SuffixArray(std::string_view s) {
    std::vector<int> text(s.size() + 1);
    std::copy(s.cbegin(), s.cend(), text.begin());

    {
      const int least =
          *std::min_element(text.cbegin(), text.cbegin() + s.size());
      text.back() = least - 1;
      for (int& c : text) {
        c -= least - 1;
      }
    }

    const int n = static_cast<int>(text.size());
    const int alpha = 1 + *std::max_element(text.cbegin(), text.cend());
    const int buckets = std::max(n, alpha) + 1;

    sa_.resize(n);
    rank_.resize(n);
    std::vector<int> aux(n);
    std::vector<int> count(buckets);

    // aux - current order
    // sa - destination
    auto CountingSort = [&](const std::vector<int>& key) -> void {
      for (int i : aux) {
        count[key[i] + 1] += 1;
      }
      for (int i = 1; i < buckets; ++i) {
        count[i] += count[i - 1];
      }
      for (int i : aux) {
        sa_[count[key[i]]++] = i;
      }
      std::fill(count.begin(), count.end(), 0);
    };

    std::iota(aux.begin(), aux.end(), 0);
    CountingSort(text);

    int current_rank = 0;
    for (int i = 1; i < n; ++i) {
      current_rank += (text[sa_[i - 1]] != text[sa_[i]]);
      rank_[sa_[i]] = current_rank;
    }

    int length = 1;
    while (current_rank != n - 1) {
      for (int i = 0; i < n; ++i) {
        if ((aux[i] = sa_[i] - length) < 0) {
          aux[i] += n;
        }
      }
      CountingSort(rank_);

      current_rank = 0;
      aux[sa_.front()] = 0;
      for (int i = 1; i < n; ++i) {
        if (rank_[sa_[i - 1]] == rank_[sa_[i]]) {
          int x = sa_[i - 1] + length < n ? rank_[sa_[i - 1] + length] : -1;
          int y = sa_[i] + length < n ? rank_[sa_[i] + length] : -1;
          current_rank += (x != y);
        } else {
          current_rank += 1;
        }
        aux[sa_[i]] = current_rank;
      }

      length *= 2;
      std::swap(aux, rank_);
    }

    sa_.erase(sa_.begin());
    rank_.pop_back();
    for (int& r : rank_) {
      r -= 1;
    }

    BuildLCP(s);
  }

  void BuildLCP(std::string_view s) {
    const int n = static_cast<int>(sa_.size());
    lcp_.resize(n);
    for (int i = 0, k = 0; i < n; ++i) {
      if (k > 0) {
        k -= 1;
      }
      if (rank_[i] == n - 1) {
        k = 0;
        continue;
      }
      int j = sa_[rank_[i] + 1];
      while (i + k < n && j + k < n && s[i + k] == s[j + k]) {
        k += 1;
      }
      lcp_[rank_[i]] = k;
    }
  }

  [[nodiscard]] const std::vector<int>& sa() const { return sa_; }
  [[nodiscard]] const std::vector<int>& rank() const { return rank_; }
  [[nodiscard]] const std::vector<int>& lcp() const { return lcp_; }

 private:
  std::vector<int> sa_;
  std::vector<int> lcp_;
  std::vector<int> rank_;
};
