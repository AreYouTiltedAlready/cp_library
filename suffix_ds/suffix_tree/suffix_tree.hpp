#include <string>
#include <vector>

namespace suffix_ds {
namespace suffix_tree {
class SuffixTree {
 public:
  struct Node {
    Node() : s_start(0), length(0), parent(0), weighted_depth(0) {}

    std::vector<std::pair<int, int>> edges;
    int parent;
    int s_start;
    int length;
    int weighted_depth;
  };

  explicit SuffixTree(std::string s, const std::vector<int>& sa,
                      const std::vector<int>& lcp)
      : s_(std::move(s)),
        tree_(s.length() * 2 + 2),
        n_(static_cast<int>(s.length())) {
    int last = MakeNode(1, sa[0], n_ - sa[0]);
    tree_[1].edges.emplace_back(s_[sa[0]], last);
    for (int i = 1; i < n_; ++i) {
      last = Extend(last, sa[i], lcp[i - 1]);
    }
    int id = 0;
    int euler_id = 0;
    tour_list_.resize(size_ - 1);
    euler_tour_.resize(2 * size_ - 1);
    auto Dfs = [&](auto& self, int v) -> void {
      tour_list_[id++] = v;
      euler_tour_[euler_id++] = v;
      for (const auto& [c, to] : tree_[v].edges) {
        self(self, to);
        euler_tour_[euler_id++] = v;
      }
    };
    Dfs(Dfs, 1);
  }

 private:
  int MakeNode(int parent, int s_start, int length) {
    Node& node = tree_[size_++];
    node.parent = parent;
    node.s_start = s_start;
    node.length = length;
    node.weighted_depth = tree_[parent].weighted_depth + length;
    return size_ - 1;
  }

  int Split(int state, int e_pos) {
    Node& current = tree_[state];
    Node& parent = tree_[current.parent];

    if (e_pos == current.length) {
      return state;
    }

    int s_start = current.s_start;
    current.length -= e_pos;
    current.s_start += e_pos;
    int result = MakeNode(current.parent, s_start, e_pos);
    parent.edges.back().second = current.parent = result;
    tree_[result].edges.emplace_back(s_[s_start + e_pos], state);
    return result;
  }

  int Extend(int state, int id, int lcp) {
    int weighted_depth = tree_[state].weighted_depth;
    int e_pos = tree_[state].length;
    while (weighted_depth > lcp) {
      if (int delta = weighted_depth - lcp; e_pos <= delta) {
        weighted_depth -= e_pos;
        state = tree_[state].parent;
        e_pos = tree_[state].length;
      } else {
        weighted_depth -= delta;
        e_pos -= delta;
      }
    }
    id += lcp;
    state = Split(state, e_pos);
    int next = MakeNode(state, id, n_ - id);
    tree_[state].edges.emplace_back(s_[id], next);
    return next;
  }

  std::vector<Node> tree_;
  std::vector<int> tour_list_;
  std::vector<int> euler_tour_;
  std::string s_;

  int size_{};
  int n_;
};

}  // namespace suffix_tree
}  // namespace suffix_ds
