#include <algorithm>
#include <optional>
#include <vector>

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

template <typename T>
concept undirected_graph =
    has_edge_t<T> &&
    std::is_base_of_v<GraphBase<typename T::edge_t, GraphType::kUndirected>, T>;

template <typename T>
concept directed_graph =
    has_edge_t<T> &&
    std::is_base_of_v<GraphBase<typename T::edge_t, GraphType::kDirected>, T>;

template <typename T>
concept weighted_graph = has_edge_t<T> && weighted_edge<typename T::edge_t>;

template <typename T>
concept graph = undirected_graph<T> || directed_graph<T>;

}  // namespace internal

using Undigraph = internal::Graph<internal::BasicEdge>;
using Digraph = internal::Digraph<internal::BasicEdge>;

template <typename Cost>
using WeightedUndigraph = internal::Graph<internal::WeightedEdge<Cost>>;

template <typename Cost>
using WeightedDigraph = internal::Digraph<internal::WeightedEdge<Cost>>;

}  // namespace graph_traits

}  // namespace graphs
