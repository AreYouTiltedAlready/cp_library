#include <algorithm>
#include <array>
#include <numeric>
#include <string_view>
#include <vector>

std::array<std::vector<int>, 2> SuffixArray(std::string_view s) {
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

  std::vector<int> sa(n);
  std::vector<int> aux(n);
  std::vector<int> rank(n);
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
      sa[count[key[i]]++] = i;
    }
    std::fill(count.begin(), count.end(), 0);
  };

  std::iota(aux.begin(), aux.end(), 0);
  CountingSort(text);

  int current_rank = 0;
  for (int i = 1; i < n; ++i) {
    if (text[sa[i]] != text[sa[i - 1]]) {
      current_rank += 1;
    }
    rank[sa[i]] = current_rank;
  }

  int length = 1;
  while (current_rank != n - 1) {
    for (int i = 0; i < n; ++i) {
      if ((aux[i] = sa[i] - length) < 0) {
        aux[i] += n;
      }
    }

    CountingSort(rank);
    current_rank = 0;
    aux[sa.front()] = 0;
    for (int i = 1; i < n; ++i) {
      int l = sa[i - 1];
      int r = sa[i];
      if (rank[l] == rank[r]) {
        if ((l += length) >= n) {
          l -= n;
        }
        if ((r += length) >= n) {
          r -= n;
        }
        current_rank += rank[l] != rank[r];
      } else {
        current_rank += 1;
      }
      aux[sa[i]] = current_rank;
    }

    length *= 2;
    std::swap(aux, rank);
  }

  sa.erase(sa.begin());
  rank.pop_back();
  for (int& r : rank) {
    r -= 1;
  }

  return std::array{sa, rank};
}

std::vector<int> LCP(std::string_view s, const std::vector<int>& sa,
                     const std::vector<int>& rank) {
  const int n = static_cast<int>(s.size());
  std::vector<int> lcp(n);
  for (int i = 0, k = 0; i < n; ++i) {
    if (k > 0) {
      k -= 1;
    }
    if (rank[i] == n - 1) {
      k = 0;
      continue;
    }
    int j = sa[rank[i] + 1];
    while (i + k < n && j + k < n && s[i + k] == s[j + k]) {
      k += 1;
    }
    lcp[rank[i]] = k;
  }
  return lcp;
}
