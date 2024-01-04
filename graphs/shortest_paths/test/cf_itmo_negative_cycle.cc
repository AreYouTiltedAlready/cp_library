// Problem: https://codeforces.com/group/dAhOSPf3oD/contest/408564/problem/C
// Submission:
// https://codeforces.com/group/dAhOSPf3oD/contest/408564/submission/240164279

#include <chrono>
#include <ext/pb_ds/priority_queue.hpp>
#include <iostream>
#include <numeric>
#include <random>

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

namespace utils {

template <typename T>
class simple_queue {
 public:
  simple_queue() : pos_(0), payload_({}) {}

  explicit simple_queue(int n) : simple_queue() { payload_.reserve(n); }

  void reserve(int n) { payload_.reserve(n); }

  void push(const T& value) { payload_.push_back(value); }

  void push(T&& value) { payload_.push_back(std::forward<T>(value)); }

  T poll() { return payload_[pos_++]; }

  [[nodiscard]] bool empty() const {
    return pos_ == static_cast<int>(payload_.size());
  }

  [[nodiscard]] int size() const {
    return static_cast<int>(payload_.size()) - pos_;
  }

 private:
  int pos_;
  std::vector<T> payload_;
};

}  // namespace utils

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

  [[nodiscard]] const Edge& edge(int id) const { return edges_[id]; }

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

namespace shortest_paths {

template <typename Cost>
struct SSSPInfo {
  using cost_t = Cost;
  static constexpr cost_t kUnreachableSentinel =
      std::numeric_limits<cost_t>::max() / 2;
  explicit SSSPInfo(int n) : distance(n, kUnreachableSentinel), p_edge(n, -1) {}

  std::vector<Cost> distance;
  std::vector<int> p_edge;
};

template <typename Cost>
struct FordBellmanInfo {
  explicit FordBellmanInfo(int n) : sp_info(n) {}

  SSSPInfo<Cost> sp_info;
  std::vector<int> negative_cycle;
};

template <typename Distance, graph_traits::weighted_graph Graph>
SSSPInfo<Distance> Dijkstra(const Graph& graph,
                            const std::vector<int>& source) {
  using distance_t = Distance;
  static constexpr distance_t kUnreachableSentinel =
      std::numeric_limits<distance_t>::max() / 2;

  SSSPInfo<distance_t> result(graph.n());
  __gnu_pbds::priority_queue<std::pair<distance_t, int>, std::greater<>> heap;
  std::vector<typename decltype(heap)::point_iterator> its(graph.n(),
                                                           heap.end());

  for (int s : source) {
    result.distance[s] = 0;
    its[s] = heap.push(std::make_pair(0, s));
  }

  const auto& g = graph.g();
  const auto& edges = graph.edges();
  std::vector<distance_t>& dist = result.distance;
  std::vector<int>& p_edge = result.p_edge;

  while (!heap.empty()) {
    auto [d, v] = heap.top();
    its[v] = heap.end();
    heap.pop();
    for (int e_id : g[v]) {
      const auto& edge = edges[e_id];
      int to = v ^ edge.from ^ edge.to;
      if (distance_t new_cost = d + edge.cost; new_cost < dist[to]) {
        p_edge[to] = e_id;
        dist[to] = new_cost;
        auto new_state = std::make_pair(new_cost, to);
        if (its[to] == heap.end()) {
          its[to] = heap.push(new_state);
        } else {
          heap.modify(its[to], new_state);
        }
      }
    }
  }

  return result;
}

template <typename Distance, graph_traits::weighted_graph Graph>
  requires graph_traits::directed_graph<Graph>
SSSPInfo<Distance> Dijkstra(const Graph& graph, int source) {
  return Dijkstra<Distance, Graph>(graph, std::vector<int>(1, source));
}

template <typename Distance, graph_traits::weighted_graph Graph>
FordBellmanInfo<Distance> FordBellman(const Graph& graph,
                                      const std::vector<int>& source) {
  using distance_t = Distance;
  FordBellmanInfo<distance_t> result(graph.n());

  for (int s : source) {
    result.sp_info.distance[s] = 0;
  }

  auto& sp_info = result.sp_info;
  auto Relax = [&graph = std::as_const(graph), &sp_info](int e_id) -> bool {
    const auto& edge = graph.edge(e_id);
    if (distance_t new_cost = sp_info.distance[edge.from] + edge.cost;
        new_cost < sp_info.distance[edge.to]) {
      sp_info.p_edge[edge.to] = e_id;
      sp_info.distance[edge.to] = new_cost;
      return true;
    }
    return false;
  };

  static std::mt19937_64 rng(
      std::chrono::steady_clock::now().time_since_epoch().count());
  std::vector<int> order(graph.m());
  std::iota(order.begin(), order.end(), 0);
  std::shuffle(order.begin(), order.end(), rng);

  bool any_relaxation = true;
  for (int i = 0; i < graph.n() - 1 && any_relaxation; ++i) {
    any_relaxation = false;
    for (int id : order) {
      any_relaxation |= Relax(id);
    }
  }

  if (!any_relaxation) {
    return result;
  }

  int last_edge = [&order = std::as_const(order), &Relax]() -> int {
    for (int id : order) {
      if (Relax(id)) {
        return id;
      }
    }
    return -1;
  }();

  if (last_edge == -1) {
    return result;
  }

  std::vector<char> seen(graph.m());
  seen[last_edge] = true;

  {
    std::vector<int> edges_stack;
    edges_stack.reserve(2 * graph.m());
    edges_stack.push_back(last_edge);
    int e_id = sp_info.p_edge[graph.edge(last_edge).from];
    while (!seen[e_id]) {
      seen[e_id] = true;
      edges_stack.push_back(e_id);
      e_id = sp_info.p_edge[graph.edge(e_id).from];
    }
    while (edges_stack.back() != e_id) {
      result.negative_cycle.push_back(edges_stack.back());
      edges_stack.pop_back();
    }
    result.negative_cycle.push_back(e_id);
  }

  return result;
}

template <typename Distance, graph_traits::weighted_graph Graph>
FordBellmanInfo<Distance> FordBellman(const Graph& graph, int source) {
  return FordBellman<Distance, Graph>(graph, source);
}

template <graph_traits::graph Graph>
SSSPInfo<int> Bfs(const Graph& graph, const std::vector<int>& source) {
  SSSPInfo<int> result(graph.n());
  ::utils::simple_queue<int> q;
  for (int s : source) {
    q.push(s);
    result.distance[s] = 0;
  }

  const auto& g = graph.g();
  const auto& edges = graph.edges();
  std::vector<int>& dist = result.distance;
  std::vector<int>& p_edge = result.p_edge;

  while (!q.empty()) {
    int v = q.poll();
    for (int e_id : g[v]) {
      const auto& edge = edges[e_id];
      int to = v ^ edge.from ^ edge.to;
      if (dist[to] == SSSPInfo<int>::kUnreachableSentinel) {
        q.push(to);
        p_edge[to] = e_id;
        dist[to] = dist[v] + 1;
      }
    }
  }

  return result;
}

template <graph_traits::graph Graph>
SSSPInfo<int> Bfs(const Graph& graph, int source) {
  return Bfs(graph, std::vector<int>(1, source));
}

}  // namespace shortest_paths

}  // namespace graphs

void RunCase([[maybe_unused]] int testcase) {
  int n;
  std::cin >> n;

  graphs::WeightedDigraph<int> g(n, n * n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      int cost;
      std::cin >> cost;
      if (cost != 100000) /* Ok.. */ {
        g.AddEdge(i, j, cost);
      }
    }
  }

  std::vector<int> source(n);
  std::iota(source.begin(), source.end(), 0);
  auto result = graphs::shortest_paths::FordBellman<int>(g, source);

  if (result.negative_cycle.empty()) {
    std::cout << "NO\n";
    return;
  }

  std::cout << "YES\n" << result.negative_cycle.size() << "\n";
  for (int id : result.negative_cycle) {
    std::cout << g.edge(id).from + 1 << " ";
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
