// Problem: https://codeforces.com/group/dAhOSPf3oD/contest/404821/problem/D
// Submission: https://codeforces.com/group/dAhOSPf3oD/contest/404821/submission/240489381

#include <bits/stdc++.h>

#include <limits>
#include <utility>

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

namespace graphs {

namespace graph_traits {
struct BasicEdge {
  BasicEdge() : from(0), to(0) {}
  BasicEdge(int from, int to) : from(from), to(to) {}

  int from;
  int to;
};

template <typename T>
requires std::is_arithmetic_v<T>
struct WeightedEdge : public BasicEdge {
  using cost_t = T;

  WeightedEdge() : BasicEdge(), cost(0) {}
  WeightedEdge(int from, int to, T cost) : BasicEdge(from, to), cost(cost) {}

  T cost;
};

template <typename T>
concept basic_edge = std::same_as<T, BasicEdge>;

template <typename T>
concept weighted_edge =
    std::is_base_of_v<BasicEdge, T> &&
    std::same_as < std::enable_if_t < std::is_arithmetic_v<typename T::cost_t>,
void >, void > ;

template <typename T>
concept edge = basic_edge<T> || weighted_edge<T>;

enum class GraphType { kDirected, kUndirected };

template <edge Edge, GraphType kGraphType>
class GraphBase {
 public:
  [[nodiscard]] int n() const { return n_; }

  [[nodiscard]] int m() const { return m_; }

  [[nodiscard]] const Edge& edge(int id) const { return edges_[id]; }

  Edge* mutable_edge(int id) { return &(edges_[id]); }

  [[nodiscard]] const std::vector<Edge>& edges() const { return edges_; }

  [[nodiscard]] const std::vector<std::vector<int>>& g() const { return g_; }

 protected:
  explicit GraphBase(int n, int m) : g_(n), edges_(), n_(n), m_(0) {
    edges_.reserve(m);
  }

  std::vector<Edge> edges_;
  std::vector<std::vector<int>> g_;

  const int n_;
  int m_;
};

template <edge Edge>
class Digraph : public GraphBase<Edge, GraphType::kDirected> {
  using base = GraphBase<Edge, GraphType::kDirected>;
  using base::edges_;
  using base::g_;
  using base::m_;
  using base::n_;

 public:
  using edge_t = Edge;

  explicit Digraph(int n, int m) : base(n, m) {}

  template <typename... Args>
  requires std::is_constructible_v<edge_t, Args...>
  void AddEdge(Args&&... args) {
    int id = m_++;
    edges_.emplace_back(std::forward<Args>(args)...);
    g_[edges_.back().from].push_back(id);
  }
};

template <edge Edge>
class Graph : public GraphBase<Edge, GraphType::kUndirected> {
  using base = GraphBase<Edge, GraphType::kUndirected>;
  using base::edges_;
  using base::g_;
  using base::m_;
  using base::n_;

 public:
  using edge_t = Edge;
  explicit Graph(int n, int m) : base(n, m) {}

  template <typename... Args>
  requires std::is_constructible_v<edge_t, Args...>
  void AddEdge(Args&&... args) {
    int id = m_++;
    edges_.emplace_back(std::forward<Args>(args)...);
    g_[edges_.back().from].push_back(id);
    g_[edges_.back().to].push_back(id);
  }
};

template <typename T>
concept has_edge_t = std::same_as < std::enable_if_t < edge<typename T::edge_t>,
void >, void > ;

template <typename T>
concept undirected_graph = has_edge_t<T> &&
    std::is_base_of_v<GraphBase<typename T::edge_t, GraphType::kUndirected>, T>;

template <typename T>
concept directed_graph = has_edge_t<T> &&
    std::is_base_of_v<GraphBase<typename T::edge_t, GraphType::kDirected>, T>;

template <typename T>
concept weighted_graph = has_edge_t<T> && weighted_edge<typename T::edge_t>;

template <typename T>
concept graph = undirected_graph<T> || directed_graph<T>;

}  // namespace graph_traits

using Undigraph = graph_traits::Graph<graph_traits::BasicEdge>;
using Digraph = graph_traits::Digraph<graph_traits::BasicEdge>;

template <typename Cost>
using WeightedUndigraph = graph_traits::Graph<graph_traits::WeightedEdge<Cost>>;

template <typename Cost>
using WeightedDigraph = graph_traits::Digraph<graph_traits::WeightedEdge<Cost>>;

namespace dfs_forest {
template <graph_traits::undirected_graph Graph>
class DfsForest {
 public:
  template <typename U>
  requires std::is_same_v<std::decay_t<U>, Graph>
  explicit DfsForest(U&& graph)
      : graph_(std::forward<U>(graph)),
        low_(graph_.n()),
        depth_(graph_.n(), -1),
        par_edge_(graph_.n(), -1) {
    Build();
  }

  [[nodiscard]] bool IsBridge(int e_id) const {
    const auto& edge = graph_.edge(e_id);
    auto [u, v] = std::make_pair(edge.from, edge.to);
    if (depth_[u] > depth_[v]) {
      std::swap(u, v);
    }
    return par_edge_[v] == e_id && low_[v] > depth_[u];
  }

  [[nodiscard]] bool IsArticulationPoint(int v) const {
    int children = 0;
    bool dangling_child = false;
    for (int e_id : graph_.g()[v]) {
      const auto& edge = graph_.edge(e_id);
      int to = v ^ edge.from ^ edge.to;
      if (par_edge_[to] == e_id) {
        children += 1;
        dangling_child |= low_[to] >= depth_[v];
      }
    }
    if (par_edge_[v] == -1) {
      return children > 1;
    }
    return dangling_child;
  }

  std::vector<int> EdgeBiconnectivityComponents() {
    int last_id = 0;
    std::vector<int> id(graph_.n(), -1);
    auto Dfs = [&](auto& self, int v) -> void {
      for (int e_id : graph_.g()[v]) {
        const auto& edge = graph_.edge(e_id);
        int to = v ^ edge.from ^ edge.to;
        if (id[to] != -1) {
          continue;
        }
        if (low_[to] > depth_[v]) {
          id[to] = last_id++;
        } else {
          id[to] = id[v];
        }
        self(self, to);
      }
    };

    for (int i = 0; i < graph_.n(); ++i) {
      if (id[i] == -1) {
        id[i] = last_id++;
        Dfs(Dfs, i);
      }
    }

    return id;
  }

  [[nodiscard]] std::vector<int> VertexBiconnectivityComponents() const {
    int last_id = 0;
    std::vector<int> id(graph_.m(), -1);
    std::vector<bool> used(graph_.n());
    auto Dfs = [&](auto& self, int v, int p_id) -> void {
      used[v] = true;
      for (int e_id : graph_.g()[v]) {
        if (id[e_id] != -1) {
          continue;
        }
        const auto& edge = graph_.edge(e_id);
        int to = v ^ edge.from ^ edge.to;
        if (!used[to]) {
          if (low_[to] >= depth_[v]) {
            id[e_id] = last_id++;
          } else {
            id[e_id] = id[p_id];
          }
          self(self, to, e_id);
        } else {
          id[e_id] = id[p_id];
        }
      }
    };
    for (int i = 0; i < graph_.n(); ++i) {
      if (par_edge_[i] == -1) {
        Dfs(Dfs, i, -1);
      }
    }
    return id;
  }

 private:
  void Build() {
    for (int i = 0; i < graph_.n(); ++i) {
      if (depth_[i] == -1) {
        depth_[i] = 0;
        Dfs(i);
      }
    }
  }

  void Dfs(int v) {
    low_[v] = depth_[v];
    for (int e_id : graph_.g()[v]) {
      if (e_id == par_edge_[v]) {
        continue;
      }
      const auto& edge = graph_.edge(e_id);
      int to = v ^ edge.from ^ edge.to;
      if (depth_[to] == -1) {
        par_edge_[to] = e_id;
        depth_[to] = depth_[v] + 1;
        Dfs(to);
        low_[v] = std::min(low_[v], low_[to]);
      } else {
        low_[v] = std::min(low_[v], depth_[to]);
      }
    }
  }

  Graph graph_;
  std::vector<int> low_;
  std::vector<int> depth_;
  std::vector<int> par_edge_;
};

}  // namespace dfs_forest

}  // namespace graphs

void RunCase([[maybe_unused]] int testcase) {
  int n;
  int m;
  std::cin >> n >> m;

  graphs::Undigraph graph(n, m);
  for (int i = 0; i < m; ++i) {
    int u;
    int v;
    std::cin >> u >> v;
    --u;
    --v;
    graph.AddEdge(u, v);
  }

  graphs::dfs_forest::DfsForest<graphs::Undigraph> dfs_forest(std::move(graph));
  std::vector<int> id = dfs_forest.EdgeBiconnectivityComponents();
  std::cout << 1 + *std::max_element(id.cbegin(), id.cend()) << "\n";
  for (int i : id) {
    std::cout << i + 1 << " ";
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
