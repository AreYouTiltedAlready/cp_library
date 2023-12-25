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

#include <algorithm>
#include <cassert>
#include <numeric>
#include <string>
#include <vector>

// Thanks to AtCoder library - I'm too lazy to implement SA-IS by myself
// This is just AtCoder library implementation with some minor changes (which makes it consistent with my codestyle)
namespace internal {
std::vector<int> SANaive(const std::vector<int>& s) {
  const int n = static_cast<int>(s.size());
  std::vector<int> sa(n);
  std::iota(sa.begin(), sa.end(), 0);
  std::sort(sa.begin(), sa.end(), [&](int l, int r) {
    if (l == r) return false;
    while (l < n && r < n) {
      if (s[l] != s[r]) { return s[l] < s[r]; }
      l += 1;
      r += 1;
    }
    return l == n;
  });
  return sa;
}

std::vector<int> SADoubling(const std::vector<int>& s) {
  const int n = static_cast<int>(s.size());
  std::vector<int> sa(n);
  std::vector<int> aux(n);
  std::vector<int> rank = s;
  std::iota(sa.begin(), sa.end(), 0);
  for (int k = 1; k < n; k *= 2) {
    auto Compare = [&](int x, int y) {
      if (rank[x] != rank[y]) { return rank[x] < rank[y]; }
      int rx = x + k < n ? rank[x + k] : -1;
      int ry = y + k < n ? rank[y + k] : -1;
      return rx < ry;
    };
    std::sort(sa.begin(), sa.end(), Compare);
    aux[sa[0]] = 0;
    for (int i = 1; i < n; i++) {
      aux[sa[i]] = aux[sa[i - 1]] + (Compare(sa[i - 1], sa[i]) ? 1 : 0);
    }
    std::swap(aux, rank);
  }
  return sa;
}

// SA-IS, linear-time suffix array construction
// Reference:
// G. Nong, S. Zhang, and W. H. Chan,
// Two Efficient Algorithms for Linear Time Suffix Array Construction
template <int THRESHOLD_NAIVE = 10, int THRESHOLD_DOUBLING = 40>
std::vector<int> SA_IS(const std::vector<int>& s, int upper) {
  const int n = static_cast<int>(s.size());
  if (n == 0) return {};
  if (n == 1) return {0};
  if (n == 2) {
    if (s[0] < s[1]) {
      return {0, 1};
    } else {
      return {1, 0};
    }
  }
  if (n < THRESHOLD_NAIVE) { return SANaive(s); }
  if (n < THRESHOLD_DOUBLING) { return SADoubling(s); }

  std::vector<int> sa(n);
  std::vector<bool> ls(n);
  for (int i = n - 2; i >= 0; i--) {
    ls[i] = (s[i] == s[i + 1]) ? ls[i + 1] : (s[i] < s[i + 1]);
  }
  std::vector<int> sum_l(upper + 1);
  std::vector<int> sum_s(upper + 1);
  for (int i = 0; i < n; i++) {
    if (!ls[i]) {
      sum_s[s[i]] += 1;
    } else {
      sum_l[s[i] + 1] += 1;
    }
  }
  for (int i = 0; i <= upper; i++) {
    sum_s[i] += sum_l[i];
    if (i < upper) { sum_l[i + 1] += sum_s[i]; }
  }

  auto Induce = [&](const std::vector<int>& lms) {
    std::fill(sa.begin(), sa.end(), -1);
    std::vector<int> buf(upper + 1);
    std::copy(sum_s.begin(), sum_s.end(), buf.begin());
    for (int d : lms) {
      if (d == n) { continue; }
      sa[buf[s[d]]++] = d;
    }
    std::copy(sum_l.begin(), sum_l.end(), buf.begin());
    sa[buf[s[n - 1]]++] = n - 1;
    for (int i = 0; i < n; i++) {
      int v = sa[i];
      if (v >= 1 && !ls[v - 1]) { sa[buf[s[v - 1]]++] = v - 1; }
    }
    std::copy(sum_l.begin(), sum_l.end(), buf.begin());
    for (int i = n - 1; i >= 0; i--) {
      int v = sa[i];
      if (v >= 1 && ls[v - 1]) { sa[--buf[s[v - 1] + 1]] = v - 1; }
    }
  };

  std::vector<int> lms_map(n + 1, -1);
  int m = 0;
  for (int i = 1; i < n; i++) {
    if (!ls[i - 1] && ls[i]) { lms_map[i] = m++; }
  }

  std::vector<int> lms;
  lms.reserve(m);
  for (int i = 1; i < n; i++) {
    if (!ls[i - 1] && ls[i]) { lms.push_back(i); }
  }

  Induce(lms);

  if (m) {
    std::vector<int> sorted_lms;
    sorted_lms.reserve(m);
    for (int v : sa) {
      if (lms_map[v] != -1) sorted_lms.push_back(v);
    }
    std::vector<int> rec_s(m);
    int rec_upper = 0;
    rec_s[lms_map[sorted_lms[0]]] = 0;
    for (int i = 1; i < m; i++) {
      int l = sorted_lms[i - 1];
      int r = sorted_lms[i];
      int end_l = (lms_map[l] + 1 < m) ? lms[lms_map[l] + 1] : n;
      int end_r = (lms_map[r] + 1 < m) ? lms[lms_map[r] + 1] : n;
      bool same = true;
      if (end_l - l != end_r - r) {
        same = false;
      } else {
        while (l < end_l) {
          if (s[l] != s[r]) { break; }
          l += 1;
          r += 1;
        }
        same &= (l != n && s[l] == s[r]);
      }
      if (!same) { rec_upper += 1; }
      rec_s[lms_map[sorted_lms[i]]] = rec_upper;
    }

    auto rec_sa = SA_IS<THRESHOLD_NAIVE, THRESHOLD_DOUBLING>(rec_s, rec_upper);

    for (int i = 0; i < m; i++) { sorted_lms[i] = lms[rec_sa[i]]; }
    Induce(sorted_lms);
  }
  return sa;
}

}  // namespace internal

std::vector<int> SuffixArray(const std::vector<int>& s, int upper) {
  assert(0 <= upper);
  for (int d : s) { assert(0 <= d && d <= upper); }
  return internal::SA_IS(s, upper);
}

template <typename T>
std::vector<int> SuffixArray(const std::vector<T>& s) {
  const int n = static_cast<int>(s.size());
  std::vector<int> id(n);
  iota(id.begin(), id.end(), 0);
  sort(id.begin(), id.end(), [&](int l, int r) { return s[l] < s[r]; });
  std::vector<int> s2(n);
  int now = 0;
  for (int i = 0; i < n; i++) {
    if (i > 0 && s[id[i - 1]] != s[id[i]]) { now += 1; }
    s2[id[i]] = now;
  }
  return internal::SA_IS(s2, now);
}

std::vector<int> SuffixArray(const std::string& s) {
  const int n = static_cast<int>(s.size());
  std::vector<int> s2(n);
  for (int i = 0; i < n; i++) { s2[i] = s[i]; }
  return internal::SA_IS(s2, 255);
}

// Reference:
// T. Kasai, G. Lee, H. Arimura, S. Arikawa, and K. Park,
// Linear-Time Longest-Common-Prefix Computation in Suffix Arrays and Its
// Applications
template <typename T>
std::vector<int> LCPArray(const std::vector<T>& s, const std::vector<int>& sa) {
  const int n = static_cast<int>(s.size());
  assert(n >= 1);
  std::vector<int> rank(n);
  for (int i = 0; i < n; i++) { rank[sa[i]] = i; }
  std::vector<int> lcp(n - 1);
  for (int i = 0, h = 0; i < n; i++) {
    if (h > 0) { h -= 1; }
    if (rank[i] == 0) continue;
    int j = sa[rank[i] - 1];
    while (i + h < n && j + h < n && s[i + h] == s[j + h]) { h += 1; }
    lcp[rank[i] - 1] = h;
  }
  return lcp;
}

std::vector<int> LCPArray(const std::string& s, const std::vector<int>& sa) {
  const int n = static_cast<int>(s.size());
  std::vector<int> s2(n);
  for (int i = 0; i < n; i++) { s2[i] = s[i]; }
  return LCPArray(s2, sa);
}

void RunCase([[maybe_unused]] int testcase) {
  std::string s;
  std::cin >> s;

  const std::vector<int> sa = SuffixArray(s);
  for (int id : sa) { std::cout << id << " "; }
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
