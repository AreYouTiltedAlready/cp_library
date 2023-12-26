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

// https://judge.yosupo.jp/problem/two_sat
// https://judge.yosupo.jp/submission/179608
class TwoSat {
 public:
  explicit TwoSat(int n) : n_(n), g_(2 * n) {}

  // u > 0 -> u, u < 0 -> ~u
  // v > 0 -> v, v < 0 -> ~v
  void AddClause(int u, int value_u, int v, int value_v) {
    g_[u * 2 + value_u].push_back(v * 2 + value_v);
    g_[v * 2 + (value_v ^ 1)].push_back(u * 2 + (value_u ^ 1));
  }

  void AddOrCondition(int u, int value_u, int v, int value_v) {
    AddClause(u, value_u ^ 1, v, value_v);
    AddClause(v, value_v ^ 1, u, value_u);
  }

  [[nodiscard]] std::vector<char> Solve() const {
    std::vector<char> result(n_);
    const std::vector<int> scc = FindSCC(g_);
    for (int i = 0; i < n_; ++i) {
      if (scc[i * 2] == scc[i * 2 + 1]) {
        return {};
      }
      result[i] = static_cast<char>(scc[i * 2] < scc[i * 2 + 1]);
    }
    return result;
  }

 private:
  int n_;
  std::vector<std::vector<int>> g_;
};

void RunCase([[maybe_unused]] int testcase) {
  {
    std::string s;
    std::cin >> s >> s;
  }

  int n;
  int m;
  std::cin >> n >> m;
  TwoSat ts(n);
  for (int i = 0; i < m; ++i) {
    int u;
    int v;
    int z;
    std::cin >> u >> v >> z;

    int value_u{};
    int value_v{};
    if (u < 0) {
      u = -(u + 1);
      value_u = 0;
    } else {
      u = u - 1;
      value_u = 1;
    }
    if (v < 0) {
      v = -(v + 1);
      value_v = 0;
    } else {
      v = v - 1;
      value_v = 1;
    }

    ts.AddOrCondition(u, value_u, v, value_v);
  }

  const std::vector<char> result = ts.Solve();
  if (result.empty()) {
    std::cout << "s UNSATISFIABLE\n";
  } else {
    std::cout << "s SATISFIABLE\nv ";
    for (int i = 0; i < n; ++i) {
      std::cout << (result[i] ? i + 1 : -(i + 1)) << " ";
    }
    std::cout << "0\n";
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
