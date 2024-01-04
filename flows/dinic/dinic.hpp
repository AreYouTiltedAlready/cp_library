#include <bits/stdc++.h>

namespace internal {

// simple queue on vector
// useful for constant speedup of bfs-based algorithms (e.g. dinic)
template <typename T>
class simple_queue {
 public:
  simple_queue() : pos_(0), payload_({}) {}
  explicit simple_queue(int n) : simple_queue() { reserve(n); }

  void pop() { pos_++; }
  void push(const T& value) { payload_.push_back(value); }
  void push(T&& value) { payload_.push_back(std::forward<T>(value)); }
  void reserve(int n) { payload_.reserve(n); }
  int poll() { return payload_[pos_++]; }
  int& front() { return payload_[pos_]; }

  [[nodiscard]] int size() const noexcept {
    return static_cast<int>(payload_.size()) - pos_;
  }

  [[nodiscard]] bool empty() const noexcept {
    return pos_ == static_cast<int>(payload_.size());
  }

 private:
  int pos_;
  std::vector<int> payload_;
};

}  // namespace internal

template <typename T>
class FlowGraph {
 public:
  using flow_t = T;

  struct Edge {
    Edge() : cap(0), flow(0), from(0), to(0) {}
    Edge(int from, int to, T cap) : cap(cap), flow(0), from(from), to(to) {}

    flow_t cap;
    flow_t flow;
    int from;
    int to;
  };

  explicit FlowGraph(int n, int m, int source, int sink)
      : edges_(m * 2),
        edge_ptr_(n),
        distance_(n),
        g_(n),
        residual_limit_(1),
        n_(n),
        m_(0),
        source_(source),
        sink_(sink) {}

  int AddEdge(int from, int to, int forward_cap, int backward_cap = 0) {
    int id = m_;
    edges_[m_++] = {from, to, forward_cap};
#pragma clang diagnostic push
#pragma ide diagnostic ignored "ArgumentSelectionDefects"
    edges_[m_++] = {to, from, backward_cap};
#pragma clang diagnostic pop
    g_[from].push_back(id);
    g_[to].push_back(id + 1);
    return id;
  }

  flow_t MaxFlowWithScaling() {
    for (const Edge& edge : edges_) {
      while (residual_limit_ < edge.cap) {
        residual_limit_ *= 2;
      }
    }
    flow_t max_flow = 0;
    while (residual_limit_ > 0) {
      max_flow += MaxFlow();
      residual_limit_ /= 2;
    }
    return max_flow;
  }

  flow_t MaxFlow(flow_t flow_limit = std::numeric_limits<flow_t>::max()) {
    flow_t max_flow = 0;
    while (Bfs() && max_flow < flow_limit) {
      flow_t augment = 0;
      do {
        augment = Dfs(source_, flow_limit - max_flow);
        max_flow += augment;
      } while (augment > 0 && max_flow < flow_limit);
    }
    return max_flow;
  }

  std::vector<std::pair<flow_t, std::vector<int>>> Decomposition() {
    int iteration = 0;
    std::vector<int> edges;
    std::vector<int> visit_time(n_, -1);
    auto Dfs = [&](auto& self, int v, flow_t flow_min) -> flow_t {
      if (v == sink_) {
        return flow_min;
      }
      visit_time[v] = iteration;
      for (int e_id : g_[v]) {
        Edge& edge = edges_[e_id];
        Edge& back_edge = edges_[e_id ^ 1];
        if (visit_time[edge.to] != iteration && edge.flow > 0) {
          if (flow_t result =
                  self(self, edge.to, std::min(edge.flow, flow_min));
              result > 0) {
            edge.flow -= result;
            back_edge.flow += result;
            edges.push_back(e_id);
            return result;
          }
        }
      }
      return 0;
    };

    flow_t flow = 0;
    std::vector<std::pair<int, std::vector<int>>> decomposition;
    while ((flow = Dfs(Dfs, source_, std::numeric_limits<flow_t>::max())) > 0) {
      iteration += 1;
      std::reverse(edges.begin(), edges.end());
      decomposition.emplace_back(flow, edges);
      edges.clear();
    }

    return decomposition;
  }

  std::vector<int> MinCut() {
    internal::simple_queue<int> que(n_);
    std::vector<char> is_reachable(n_);
    que.push(source_);
    is_reachable[source_] = true;

    while (!que.empty()) {
      int v = que.poll();
      for (int e_id : g_[v]) {
        const Edge& edge = edges_[e_id];
        if (edge.flow != edge.cap && !is_reachable[edge.to]) {
          que.push(edge.to);
          is_reachable[edge.to] = true;
        }
      }
    }

    std::vector<int> result;
    result.reserve(m_ / 2);
    for (int i = 0; i < m_ / 2; ++i) {
      auto [u, v] = std::make_pair(edges_[i * 2].from, edges_[i * 2].to);
      if (is_reachable[u] && !is_reachable[v]) {
        result.push_back(i * 2);
      } else if (!is_reachable[u] && is_reachable[v]) {
        result.push_back(i * 2 + 1);
      }
    }
    return result;
  }

  [[nodiscard]] const std::vector<Edge>& Edges() const { return edges_; }

  [[nodiscard]] const Edge& GetEdge(int id) { return edges_[id]; }

 private:
  flow_t Dfs(int v, int residual_min) {
    if (v == sink_) {
      return residual_min;
    }
    for (int& i = edge_ptr_[v]; i < static_cast<int>(g_[v].size()); ++i) {
      int id = g_[v][i];
      Edge& edge = edges_[id];
      Edge& back_edge = edges_[id ^ 1];
      if (flow_t residual = edge.cap - edge.flow;
          residual >= residual_limit_ &&
          distance_[edge.to] == distance_[v] + 1) {
        if (flow_t result = Dfs(edge.to, std::min(residual, residual_min));
            result > 0) {
          edge.flow += result;
          back_edge.flow -= result;
          return result;
        }
      }
    }
    return 0;
  }

  bool Bfs() {
    internal::simple_queue<int> que(n_);
    std::fill(edge_ptr_.begin(), edge_ptr_.end(), 0);
    std::fill(distance_.begin(), distance_.end(), -1);

    que.push(source_);
    distance_[source_] = 0;
    while (!que.empty()) {
      int v = que.poll();
      for (int e_id : g_[v]) {
        const Edge& edge = edges_[e_id];
        if (edge.cap - edge.flow >= residual_limit_ &&
            distance_[edge.to] == -1) {
          que.push(edge.to);
          distance_[edge.to] = distance_[v] + 1;
        }
      }
    }

    return distance_[sink_] != -1;
  }

  std::vector<Edge> edges_;
  std::vector<int> edge_ptr_;
  std::vector<int> distance_;
  std::vector<std::vector<int>> g_;

  flow_t residual_limit_;
  int n_;
  int m_;
  const int source_;
  const int sink_;
};
