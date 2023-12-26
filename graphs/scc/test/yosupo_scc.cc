#include <bits/stdc++.h>

namespace {

template <class Fun>
class y_combinator_result {
  Fun fun_;

 public:
  template <class T>
  explicit y_combinator_result(T&& fun)
      : fun_(std::forward<T>(fun)) {
  }  // NOLINT(*-forwarding-reference-overload)

  template <class... Args>
  decltype(auto) operator()(Args&&... args) {
    return fun_(std::ref(*this), std::forward<Args>(args)...);
  }
};

template <class Fun>
decltype(auto) y_combinator(Fun&& fun) {
  return y_combinator_result<std::decay_t<Fun>>(std::forward<Fun>(fun));
}

std::vector<int> FindSCC(const std::vector<std::vector<int>>& g) {
  const int n = static_cast<int>(g.size());
  std::vector<int> order;
  order.reserve(n);
  {
    std::vector<char> used(n);
    auto Dfs = y_combinator([&](auto&& dfs, int v) -> void {
      used[v] = true;
      for (int to : g[v]) {
        if (!used[to]) {
          dfs(to);
        }
      }
      order.push_back(v);
    });

    for (int i = 0; i < n; ++i) {
      if (!used[i]) {
        Dfs(i);
      }
    }
  }
  const auto g_rev = [](const std::vector<std::vector<int>>& g)
      -> std::vector<std::vector<int>> {
    const int n = static_cast<int>(g.size());
    std::vector<std::vector<int>> g_rev(n);
    for (int i = 0; i < n; ++i) {
      for (int to : g[i]) {
        g_rev[to].push_back(i);
      }
    }
    return g_rev;
  }(g);
  int current_id = 0;
  std::vector<int> scc_id(n, -1);
  auto Dfs = y_combinator([&](auto&& dfs, int v) -> void {
    scc_id[v] = current_id;
    for (int to : g_rev[v]) {
      if (scc_id[to] == -1) {
        dfs(to);
      }
    }
  });
  std::reverse(order.begin(), order.end());
  for (int id : order) {
    if (scc_id[id] == -1) {
      Dfs(id);
      current_id += 1;
    }
  }
  return scc_id;
}

// https://judge.yosupo.jp/problem/scc
// https://judge.yosupo.jp/submission/179604
void RunCase([[maybe_unused]] int testcase) {
  int n;
  int m;
  std::cin >> n >> m;

  std::vector<std::vector<int>> g(n);
  for (int i = 0; i < m; ++i) {
    int u;
    int v;
    std::cin >> u >> v;
    g[u].push_back(v);
  }

  std::vector<int> scc_id = FindSCC(g);
  const int k = 1 + *std::max_element(scc_id.cbegin(), scc_id.cend());
  std::vector<int> count(k);
  for (int i = 0; i < n; ++i) {
    count[scc_id[i]] += 1;
  }
  std::vector<std::vector<int>> nodes(k);
  for (int i = 0; i < k; ++i) {
    nodes[i].reserve(count[i]);
  }
  for (int i = 0; i < n; ++i) {
    nodes[scc_id[i]].push_back(i);
  }
  std::cout << k << "\n";
  for (int i = 0; i < k; ++i) {
    std::cout << count[i] << " ";
    for (int v : nodes[i]) {
      std::cout << v << " ";
    }
    std::cout << "\n";
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
