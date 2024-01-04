// Problem: https://judge.yosupo.jp/problem/scc
// Submission: https://judge.yosupo.jp/submission/181391

#include <ext/pb_ds/priority_queue.hpp>
#include <iostream>
#include <optional>

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
namespace internal {

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
    std::same_as<
        std::enable_if_t<std::is_arithmetic_v<typename T::cost_t>, void>, void>;

template <typename T>
concept edge = basic_edge<T> || weighted_edge<T>;

enum class GraphType { kDirected, kUndirected };

template <edge Edge, GraphType kGraphType>
class GraphBase {
 public:
  [[nodiscard]] int n() const { return n_; }

  [[nodiscard]] int m() const { return m_; }

  [[nodiscard]] const Edge& edge(int id) { return edges_[id]; }

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
concept has_edge_t =
    std::same_as<std::enable_if_t<edge<typename T::edge_t>, void>, void>;

}  // namespace internal

template <typename T>
concept undirected_graph =
    internal::has_edge_t<T> &&
    std::is_base_of_v<internal::GraphBase<typename T::edge_t,
                                          internal::GraphType::kUndirected>,
                      T>;

template <typename T>
concept directed_graph =
    internal::has_edge_t<T> &&
    std::is_base_of_v<
        internal::GraphBase<typename T::edge_t, internal::GraphType::kDirected>,
        T>;

template <typename T>
concept weighted_graph =
    internal::has_edge_t<T> && internal::weighted_edge<typename T::edge_t>;

template <typename T>
concept graph = undirected_graph<T> || directed_graph<T>;

}  // namespace graph_traits

using Undigraph =
    graph_traits::internal::Graph<graph_traits::internal::BasicEdge>;
using Digraph =
    graph_traits::internal::Digraph<graph_traits::internal::BasicEdge>;

template <typename Cost>
using WeightedUndigraph =
    graph_traits::internal::Graph<graph_traits::internal::WeightedEdge<Cost>>;

template <typename Cost>
using WeightedDigraph =
    graph_traits::internal::Digraph<graph_traits::internal::WeightedEdge<Cost>>;

namespace scc {
namespace internal {

template <graph_traits::directed_graph Graph>
Graph Reversed(const Graph& graph) {
  Graph result(graph.n(), graph.m());
  for (auto edge : graph.edges()) {
    std::swap(edge.from, edge.to);
    result.AddEdge(std::move(edge));
  }
  return result;
}

}  // namespace internal

template <graph_traits::directed_graph Graph>
std::vector<int> FindSCC(const Graph& graph) {
  const int n = graph.n();
  std::vector<int> order;
  order.reserve(n);

  {
    const auto& g = graph.g();
    const auto& edges = graph.edges();
    std::vector<char> used(n);
    auto Dfs = [&](auto& self, int v) -> void {
      used[v] = true;
      for (int e_id : g[v]) {
        int to = edges[e_id].to;
        if (!used[to]) {
          self(self, to);
        }
      }
      order.push_back(v);
    };

    for (int i = 0; i < n; ++i) {
      if (!used[i]) {
        Dfs(Dfs, i);
      }
    }
  }

  const Graph rev_graph = internal::Reversed(graph);
  const auto& rev_g = rev_graph.g();
  const auto& rev_edges = rev_graph.edges();

  int current_id = 0;
  std::vector<int> scc_id(n, -1);
  auto Dfs = [&](auto& self, int v) -> void {
    scc_id[v] = current_id;
    for (int e_id : rev_g[v]) {
      int to = rev_edges[e_id].to;
      if (scc_id[to] == -1) {
        self(self, to);
      }
    }
  };
  std::reverse(order.begin(), order.end());
  for (int id : order) {
    if (scc_id[id] == -1) {
      Dfs(Dfs, id);
      current_id += 1;
    }
  }
  return scc_id;
}

}  // namespace scc

}  // namespace graphs

void RunCase([[maybe_unused]] int testcase) {
  int n;
  int m;
  std::cin >> n >> m;

  graphs::Digraph g(n, m);
  for (int i = 0; i < m; ++i) {
    int u;
    int v;
    std::cin >> u >> v;
    g.AddEdge(u, v);
  }

  std::vector<int> scc_id = graphs::scc::FindSCC(g);
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
