#include <limits>
#include <optional>
#include <vector>

template <typename T>
class simple_queue {
 public:
  simple_queue() : pos_(0), payload_({}) {}
  explicit simple_queue(int n) : simple_queue() { reserve(n); }

  void reserve(int n) { payload_.reserve(n); }
  void push(const T& value) { payload_.push_back(value); }
  void push(T&& value) { payload_.push_back(std::forward<T>(value)); }
  T poll() { return payload_[pos_++]; }
  int size() const noexcept { return static_cast<int>(payload_.size()) - pos_; }
  bool empty() const noexcept {
    return pos_ == static_cast<int>(payload_.size());
  }

 private:
  int pos_;
  std::vector<T> payload_;
};

// Contains all what one needs for flow problems
// Dinic() finds max flow in $O(V^2E)$, DinicWithScaling() - in $O(VE\log C)$,
// where C is the max capacity among all edges
// Suppose we have a flow $F$:
// * FlowDecomposition() decomposes $F$ into cycles and $s \rightarrow t$ paths
// in $O(VE)$
// * MinCut() finds minimal cut between $s$ and $t$ in $O(E)$, yields indices of
// edges Note: $k$-Flow works in $O(kE)$, same as Ford-Fulkerson
template <typename T,
          typename std::enable_if_t<
              std::is_signed_v<T> && std::is_integral_v<T>, void>* = nullptr>
class FlowGraph {
 public:
  struct Edge {
    Edge() : from(0), to(0), cap(0), flow(0) {}
    Edge(int from, int to, T cap) : cap(cap), flow(0), from(from), to(to) {}

    T cap;
    T flow;
    int from;
    int to;
  };

  int AddEdge(int from, int to, T capacity, T rev_capacity = 0) {
    int id = static_cast<int>(edges_.size());
    edges_.emplace_back(from, to, capacity);
    edges_.emplace_back(to, from, rev_capacity);
    g_[from].push_back(id);
    g_[to].push_back(id + 1);
    return id;
  }

  explicit FlowGraph(int n, int m, int source, int sink)
      : n_(n), source_(source), sink_(sink), distance_(n), edge_ptr_(n), g_(n) {
    edges_.reserve(m * 2);
  }

  T FindKFlow(T k) {
    T flow = 0;
    while (flow < k && Bfs(1)) {
      while (flow < k && Dfs(source_, 1, 1) == 1) {
        flow += 1;
      }
    }
    return flow;
  }

  struct SimpleDecompositionResult {
    SimpleDecompositionResult() : vertices({}), flow(0) {}

    SimpleDecompositionResult(std::vector<int>&& vertices, T flow)
        : vertices(std::forward<std::vector<int>>(vertices)), flow(flow) {}

    std::vector<int> vertices;
    T flow;
  };

  std::vector<SimpleDecompositionResult> FlowDecomposition() {
    std::fill(edge_ptr_.begin(), edge_ptr_.end(), 0);
    std::vector<char> visited(n_);
    auto SimpleDecomposition =
        [&]() -> std::optional<SimpleDecompositionResult> {
      std::vector<int> eids;
      std::fill(visited.begin(), visited.end(), 0);
      int v = source_;
      while (!visited[v]) {
        if (v == sink_) {
          break;
        }
        for (int& i = edge_ptr_[v]; i < static_cast<int>(g_[v].size()); ++i) {
          if (const Edge& e = edges_[g_[v][i]]; e.flow > 0) {
            break;
          }
        }
        if (edge_ptr_[v] == static_cast<int>(g_[v].size())) {
          return std::nullopt;
        }
        int id = g_[v][edge_ptr_[v]];
        eids.push_back(id);

        visited[v] = true;
        v = edges_[id].to;
      }
      if (visited[v]) {
        int id = 0;
        while (edges_[eids[id]].from != v) {
          id += 1;
        }
        eids.erase(eids.begin(), eids.begin() + id);
      }
      std::vector<int> vertices;
      vertices.reserve(eids.size() + 1);
      T path_min = std::numeric_limits<T>::max();
      for (int id : eids) {
        path_min = std::min(path_min, edges_[id].flow);
        vertices.push_back(edges_[id].from);
      }
      vertices.push_back(v);
      for (int id : eids) {
        Push(id, -path_min);
      }
      return SimpleDecompositionResult(std::move(vertices), path_min);
    };

    std::vector<SimpleDecompositionResult> result;
    auto flow_part = SimpleDecomposition();
    while (flow_part.has_value()) {
      result.push_back(std::move(*flow_part));
      flow_part = SimpleDecomposition();
    }

    return result;
  }

  T DinicWithScaling() {
    T max_cap = 0;
    for (const auto& [cap, flow, from, to] : edges_) {
      max_cap = std::max(max_cap, cap);
    }
    if (max_cap == 0) {
      return 0;
    }
    T max_flow = 0;
    T bound = static_cast<T>(1) << std::__lg(max_cap);
    while (bound > 0) {
      max_flow += Dinic(bound);
      bound /= 2;
    }

    return max_flow;
  }

  T Dinic(T lower_bound = 1) {
    T max_flow = 0;
    while (Bfs(lower_bound)) {
      T augment = 0;
      while ((augment = Dfs(source_, std::numeric_limits<T>::max(),
                            lower_bound)) >= lower_bound) {
        max_flow += augment;
      }
    }
    return max_flow;
  }

  [[nodiscard]] const Edge& GetEdge(int id) const { return edges_[id]; }

  std::vector<int> MinCut() const {
    std::vector<char> reachable(n_);
    y_combinator([&](auto&& dfs, int v) -> void {
      reachable[v] = true;
      for (int eid : g_[v]) {
        if (const Edge& e = edges_[eid]; e.cap > e.flow && !reachable[e.to]) {
          dfs(e.to);
        }
      }
    })(source_);
    std::vector<int> result;
    result.reserve(edges_.size() / 2);
    for (int i = 0; i < static_cast<int>(edges_.size()); i += 2) {
      const Edge& e(edges_[i]);
      if (reachable[e.from] ^ reachable[e.to]) {
        int id = (reachable[e.from] ? i : i + 1);
        result.push_back(id);
      }
    }

    return result;
  }

 private:
  bool Bfs(T lower_bound) {
    std::fill(edge_ptr_.begin(), edge_ptr_.end(), 0);
    std::fill(distance_.begin(), distance_.end(), -1);
    simple_queue<int> que(n_);
    que.push(source_);
    distance_[source_] = 0;
    while (!que.empty()) {
      int v = que.poll();
      for (int eid : g_[v]) {
        const Edge& e = edges_[eid];
        if (e.cap - e.flow >= lower_bound && distance_[e.to] == -1) {
          que.push(e.to);
          distance_[e.to] = distance_[v] + 1;
        }
      }
    }
    return distance_[sink_] != -1;
  }

  T Dfs(int v, T least_residual, T lower_bound) {
    if (v == sink_) {
      return least_residual;
    }
    T dfs_result = 0;
    for (int& i = edge_ptr_[v]; i < static_cast<int>(g_[v].size()); ++i) {
      int eid = g_[v][i];
      if (const Edge& e = edges_[eid];
          distance_[e.to] == distance_[v] + 1 &&
          e.cap - e.flow >= lower_bound &&
          (dfs_result = Dfs(e.to, std::min(e.cap - e.flow, least_residual),
                            lower_bound)) >= lower_bound) {
        Push(eid, dfs_result);
        return dfs_result;
      }
    }
    return 0;
  }

  inline void Push(int eid, T flow) {
    edges_[eid].flow += flow;
    edges_[eid ^ 1].flow -= flow;
  }

  int n_;
  const int source_;
  const int sink_;
  std::vector<Edge> edges_;
  std::vector<int> distance_;
  std::vector<int> edge_ptr_;
  std::vector<std::vector<int>> g_;
};
