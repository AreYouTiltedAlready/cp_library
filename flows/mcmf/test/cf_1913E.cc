#include <bits/extc++.h>
#include <bits/stdc++.h>

#include <ext/pb_ds/priority_queue.hpp>
#include <functional>
#include <utility>

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

template <typename F, typename C>
class MinCostFlowGraph {
 public:
  struct Edge {
    Edge() : cap(0), flow(0), cost(0), from(0), to(0) {}
    Edge(int from, int to, F cap, C cost)
        : cap(cap), flow(0), cost(cost), from(from), to(to) {}

    F cap;
    F flow;
    C cost;
    int from;
    int to;
  };

  explicit MinCostFlowGraph(int n, int m)
      : n_(n),
        pot_(n),
        distance_(n),
        edges_(),
        sp_edge_(n),
        g_(n),
        heap_(),
        its_(n) {
    edges_.reserve(m * 2);
  }

  int AddEdge(int from, int to, F cap, C cost) {
    int id = static_cast<int>(edges_.size());
    edges_.emplace_back(from, to, cap, cost);
    edges_.emplace_back(to, from, 0, -cost);
    g_[from].push_back(id);
    g_[to].push_back(id + 1);
    return id;
  }

  [[nodiscard]] const Edge& GetEdge(int id) const { return edges_[id]; }

  void FindPath(int s, int t) {
    heap_.clear();
    std::fill(distance_.begin(), distance_.end(), kUnreachable);
    distance_[s] = 0;

    std::fill(its_.begin(), its_.end(), heap_.end());
    its_[s] = heap_.push(std::make_pair(0, s));
    while (!heap_.empty()) {
      auto [d, v] = heap_.top();
      d = pot_[v] - d;
      heap_.pop();
      for (int eid : g_[v]) {
        const Edge& e = edges_[eid];
        if (e.cap == e.flow) { continue; }
        if (C new_cost = d + e.cost - pot_[e.to]; new_cost < distance_[e.to]) {
          sp_edge_[e.to] = eid;
          distance_[e.to] = new_cost;
          auto new_pair = std::make_pair(-new_cost, e.to);
          if (its_[e.to] == heap_.end()) {
            its_[e.to] = heap_.push(new_pair);
          } else {
            heap_.modify(its_[e.to], new_pair);
          }
        }
      }
    }

    for (int i = 0; i < n_; ++i) {
      pot_[i] += distance_[i] != kUnreachable ? distance_[i] : 0;
    }
  }

  std::pair<F, C> MinCostFlow(int s, int t,
                              F flow_limit = std::numeric_limits<F>::max()) {
    F flow = 0;
    C cost = 0;

    const bool negative_edges = [&]() -> bool {
      for (int i = 0; i < edges_.size() / 2; ++i) {
        if (edges_[i * 2].cost < 0) { return true; }
      }
      return false;
    }();

    if (negative_edges) {
      std::fill(pot_.begin(), pot_.end(), kUnreachable);
      pot_[s] = 0;
      bool any = true;
      while (any) {
        any = false;
        for (int i = 0; i < static_cast<int>(edges_.size()); i += 2) {
          const Edge& e = edges_[i];
          if (pot_[e.from] == kUnreachable) { continue; }
          if (C new_cost = pot_[e.from] + e.cost; new_cost < pot_[e.to]) {
            any = true;
            pot_[e.to] = new_cost;
          }
        }
      }
      std::replace(pot_.begin(), pot_.end(), kUnreachable, static_cast<C>(0));
    }

    FindPath(s, t);
    while (flow < flow_limit && distance_[t] != kUnreachable) {
      F path_min = flow_limit - flow;
      {
        int v = t;
        while (v != s) {
          int eid = sp_edge_[v];
          const Edge& e = edges_[eid];
          path_min = std::min(path_min, e.cap - e.flow);
          v = e.from;
        }
      }
      flow += path_min;
      C new_cost = 0;
      {
        int v = t;
        while (v != s) {
          int eid = sp_edge_[v];
          Edge& e = edges_[eid];
          Edge& back = edges_[eid ^ 1];
          new_cost += e.cost;
          e.flow += path_min;
          back.flow -= path_min;
          v = e.from;
        }
      }

      cost += new_cost * path_min;
      FindPath(s, t);
    }

    return std::make_pair(flow, cost);
  }

 private:
  static constexpr C kUnreachable = std::numeric_limits<C>::max() / 2;

  const int n_;
  std::vector<C> pot_;
  std::vector<C> distance_;
  std::vector<Edge> edges_;
  std::vector<int> sp_edge_;
  std::vector<std::vector<int>> g_;
  __gnu_pbds::priority_queue<std::pair<C, int>> heap_;
  std::vector<typename decltype(heap_)::point_iterator> its_;
};

// https://codeforces.com/contest/1913/problem/E
// https://codeforces.com/contest/1913/submission/238802524
void RunCase([[maybe_unused]] int testcase) {
  int n;
  int m;
  std::cin >> n >> m;

  std::vector mat(n, std::vector<int>(m));
  for (std::vector<int>& row : mat) {
    for (int& c : row) { std::cin >> c; }
  }

  std::vector<int> row_sum(n);
  for (int& r : row_sum) { std::cin >> r; }

  std::vector<int> column_sum(m);
  for (int& c : column_sum) { std::cin >> c; }

  const int total = std::accumulate(row_sum.cbegin(), row_sum.cend(), 0);
  if (total != std::accumulate(column_sum.cbegin(), column_sum.cend(), 0)) {
    std::cout << "-1\n";
    return;
  }

  MinCostFlowGraph<int, int> Graph(n + m + 2, n * m + n + m);
  for (int i = 0; i < n; ++i) { Graph.AddEdge(n + m, i, row_sum[i], 0); }
  for (int i = 0; i < m; ++i) {
    Graph.AddEdge(i + n, n + m + 1, column_sum[i], 0);
  }

  std::vector id(n, std::vector<int>(m));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      id[i][j] = Graph.AddEdge(i, n + j, 1, 1 - mat[i][j]);
    }
  }

  auto [flow, cost] = Graph.MinCostFlow(n + m, n + m + 1);
  if (flow != total) {
    std::cout << "-1\n";
    return;
  }

  int ans = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      ans += mat[i][j] ^ Graph.GetEdge(id[i][j]).flow;
    }
  }

  std::cout << ans << "\n";
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
