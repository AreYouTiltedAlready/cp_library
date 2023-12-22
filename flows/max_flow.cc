#include <limits>
#include <optional>
#include <type_traits>
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

template <typename T,
          typename std::enable_if_t<
              std::is_signed_v<T> && std::is_integral_v<T>, void>* = nullptr>
class Network {
 public:
  struct Edge {
    Edge() : from(0), to(0), capacity(0), flow(0) {}
    Edge(int from, int to, T capacity)
        : capacity(capacity), flow(0), from(from), to(to) {}

    T Residual() const { return capacity - flow; }

    T capacity;
    T flow;
    int from;
    int to;
  };

  int AddEdge(int from, int to, T capacity, T rev_capacity = 0) {
    int id = static_cast<int>(edges_.size());
    edges_.emplace_back(from, to, capacity);
    edges_.emplace_back(to, from, rev_capacity);
    adj_list_[from].push_back(id);
    adj_list_[to].push_back(id + 1);
    return id;
  }

  explicit Network(int n, int m, int source, int sink)
      : n_(n),
        source_(source),
        sink_(sink),
        distance_(n),
        edge_ptr_(n),
        adj_list_(n) {
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
      std::vector<int> edge_ids;
      std::fill(visited.begin(), visited.end(), 0);

      int v = source_;
      while (!visited[v]) {
        if (v == sink_) {
          break;
        }

        for (int& i = edge_ptr_[v]; i < static_cast<int>(adj_list_[v].size());
             ++i) {
          if (const Edge& e = edges_[adj_list_[v][i]]; e.flow > 0) {
            break;
          }
        }

        if (edge_ptr_[v] == static_cast<int>(adj_list_[v].size())) {
          return std::nullopt;
        }
        int id = adj_list_[v][edge_ptr_[v]];
        edge_ids.push_back(id);

        visited[v] = true;
        v = edges_[id].to;
      }

      if (visited[v]) {
        int id = 0;
        while (edges_[edge_ids[id]].from != v) {
          id += 1;
        }
        edge_ids.erase(edge_ids.begin(), edge_ids.begin() + id);
      }

      std::vector<int> vertices;
      vertices.reserve(edge_ids.size() + 1);
      T path_min = std::numeric_limits<T>::max();
      for (int id : edge_ids) {
        path_min = std::min(path_min, edges_[id].flow);
        vertices.push_back(edges_[id].from);
      }

      vertices.push_back(v);
      for (int id : edge_ids) {
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
    T max_capacity = 0;
    for (const auto& [capacity, flow, from, to] : edges_) {
      max_capacity = std::max(max_capacity, capacity);
    }

    if (max_capacity == 0) {
      return 0;
    }

    T max_flow = 0;
    T bound = static_cast<T>(1) << std::__lg(max_capacity);
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

  std::vector<Edge> MinCut() const {
    std::vector<char> reachable(n_);
    y_combinator([&](auto&& dfs, int v) -> void {
      reachable[v] = true;
      for (int edge_id : adj_list_[v]) {
        if (const Edge& e = edges_[edge_id];
            e.Residual() > 0 && !reachable[e.to]) {
          dfs(e.to);
        }
      }
    })(source_);

    std::vector<Edge> result;
    result.reserve(edges_.size() / 2);
    for (int i = 0; i < static_cast<int>(edges_.size()); i += 2) {
      const Edge& e(edges_[i]);
      if (reachable[e.from] ^ reachable[e.to]) {
        int id = (reachable[e.from] ? i : i + 1);
        result.push_back(edges_[id]);
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
      for (int edge_id : adj_list_[v]) {
        const Edge& e = edges_[edge_id];
        if (e.Residual() >= lower_bound && distance_[e.to] == -1) {
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
    for (int& i = edge_ptr_[v]; i < static_cast<int>(adj_list_[v].size());
         ++i) {
      int edge_id = adj_list_[v][i];
      if (const Edge& e = edges_[edge_id];
          distance_[e.to] == distance_[v] + 1 && e.Residual() >= lower_bound &&
          (dfs_result = Dfs(e.to, std::min(e.Residual(), least_residual),
                            lower_bound)) >= lower_bound) {
        Push(edge_id, dfs_result);
        return dfs_result;
      }
    }

    return 0;
  }

  inline void Push(int edge_id, T flow) {
    edges_[edge_id].flow += flow;
    edges_[edge_id ^ 1].flow -= flow;
  }

  int n_;
  const int source_;
  const int sink_;
  std::vector<Edge> edges_;
  std::vector<int> distance_;
  std::vector<int> edge_ptr_;
  std::vector<std::vector<int>> adj_list_;
};
