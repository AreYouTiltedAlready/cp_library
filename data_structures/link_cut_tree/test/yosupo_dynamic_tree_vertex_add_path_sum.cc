// Problem: https://judge.yosupo.jp/problem/dynamic_tree_vertex_add_path_sum
// Submission: https://judge.yosupo.jp/submission/182565

#include <bits/stdc++.h>

#include <cassert>
#include <functional>

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

namespace link_cut_tree {

namespace internal {
template <typename T>
constexpr bool is_binary_addable = requires(T& lhs, const T& rhs) {
  { lhs += rhs } -> std::same_as<T&>;
};

template <typename T>
constexpr bool is_invertible = requires(T& lhs, const T& rhs) {
  { lhs -= rhs } -> std::same_as<T&>;
};

template <typename T>
concept invertible_monoid = std::is_default_constructible_v<T> &&
    is_binary_addable<T> && is_invertible<T>;

}  // namespace internal

// Overload += and -= too
// Operation must be commutative (lhs + rhs == rhs + lhs)
// Basically, you need invertibility only if you want to support subtree queries
// Otherwise, you can define - and -= (and even +=) as dummy as your imagination
// allows
template <internal::invertible_monoid S>
class LinkCutTree {
 public:
  explicit LinkCutTree(int n) : tree_(n) {}

  S GetSubtree(int parent, int child) {
    tree_[parent].MakeRoot();
    tree_[child].Access();
    return tree_[child].subtree_aggregate - tree_[parent].subtree_aggregate;
  }

  S GetPath(int u, int v) {
    tree_[u].MakeRoot();
    if (u == v) {
      return tree_[v].value;
    }
    tree_[v].Access();
    return tree_[v].left->splay_subtree_aggregate + tree_[v].value;
  }

  void Apply(int v, const S& value) {
    tree_[v].MakeRoot();
    tree_[v].value += value;
    tree_[v].Pull();
  }

  void Set(int v, const S& new_value) {
    tree_[v].MakeRoot();
    tree_[v].value = new_value;
    tree_[v].Pull();
  }

  void Link(int u, int v) { Link(&tree_[u], &tree_[v]); }

  void Cut(int u, int v) { Cut(&tree_[u], &tree_[v]); }

  bool Connected(int u, int v) { return Connected(&tree_[u], &tree_[v]); }

 private:
  struct Node {
    Node()
        : value(),
          subtree_aggregate(),
          splay_subtree_aggregate(),
          virtual_aggregate(),
          size(1) {}

    void Push() {
      if (is_reversed == 0) {
        return;
      }
      is_reversed = 0;
      std::swap(left, right);
      if (left != nullptr) {
        left->is_reversed ^= 1;
      }
      if (right != nullptr) {
        right->is_reversed ^= 1;
      }
      Pull();
    }

    void Pull() {
      size = 1;
      subtree_aggregate = virtual_aggregate + value;
      splay_subtree_aggregate = value;
      if (left != nullptr) {
        size += left->size;
        subtree_aggregate += left->subtree_aggregate;
        splay_subtree_aggregate += left->splay_subtree_aggregate;
      }
      if (right != nullptr) {
        size += right->size;
        subtree_aggregate += right->subtree_aggregate;
        splay_subtree_aggregate += right->splay_subtree_aggregate;
      }
    }

    void Rotate() {
      Node* par = this->parent;
      Node* g = par->parent;
      if (g != nullptr) {
        if (par == g->left) {
          g->left = this;
        } else {
          g->right = this;
        }
      }
      par->Push();
      this->Push();
      if (this == par->left) {
        par->left = std::exchange(this->right, par);
        if (par->left != nullptr) {
          par->left->parent = par;
        }
      } else {
        par->right = std::exchange(this->left, par);
        if (par->right != nullptr) {
          par->right->parent = par;
        }
      }

      this->parent = g;
      par->parent = this;
      par->Pull();
      this->Pull();
    }

    void Splay() {
      while (this->parent != nullptr) {
        Node* par = this->parent;
        if (par->path_parent != nullptr) {
          std::swap(this->path_parent, par->path_parent);
        }
        if (Node* g = par->parent; g != nullptr) {
          if (g->path_parent != nullptr) {
            std::swap(this->path_parent, g->path_parent);
          }
          if ((this == par->left) == (par == g->left)) {
            par->Rotate();
          } else {
            this->Rotate();
          }
        }
        this->Rotate();
      }
      this->Push();
    }

    void Access() {
      this->Splay();
      while (this->path_parent != nullptr) {
        Node* pp = this->path_parent;
        pp->Splay();
        std::swap(this->parent, this->path_parent);
        if (Node* old_right = pp->right; old_right != nullptr) {
          pp->virtual_aggregate += old_right->subtree_aggregate;
          std::swap(old_right->path_parent, old_right->parent);
        }
        pp->right = this;
        pp->virtual_aggregate -= this->subtree_aggregate;
        pp->Pull();
        this->Splay();
      }
    }

    void MakeRoot() {
      Access();
      if (right != nullptr) {
        std::swap(right->parent, right->path_parent);
        this->virtual_aggregate += right->subtree_aggregate;
        right = nullptr;
        Pull();
      }
      is_reversed ^= 1;
    }

    Node* parent{};
    Node* path_parent{};
    Node* left{};
    Node* right{};
    S value;
    S virtual_aggregate;
    S subtree_aggregate;
    S splay_subtree_aggregate;
    int size{};
    int subtree_size{};
    int is_reversed{};
  };

  static void Link(Node* u, Node* v) {
    u->MakeRoot();
    v->Access();
    u->path_parent = v;
    v->virtual_aggregate += u->subtree_aggregate;
    v->Pull();
  }

  static void Cut(Node* u, Node* v) {
    u->MakeRoot();
    v->Access();
    assert(v->left == u);
    u->parent = nullptr;
    v->left = nullptr;
    v->Pull();
  }

  static bool Connected(Node* u, Node* v) {
    u->Access();
    v->Access();
    return (u->parent != nullptr || u->path_parent != nullptr);
  }

  std::vector<Node> tree_;
};

}  // namespace link_cut_tree

void RunCase([[maybe_unused]] int testcase) {
  int n{};
  int q{};
  std::cin >> n >> q;

  link_cut_tree::LinkCutTree<int64_t> link_cut_tree(n);
  for (int i = 0; i < n; ++i) {
    int value{};
    std::cin >> value;
    link_cut_tree.Set(i, value);
  }

  for (int i = 0; i < n - 1; ++i) {
    int u{};
    int v{};
    std::cin >> u >> v;
    link_cut_tree.Link(u, v);
  }

  for (int i = 0; i < q; ++i) {
    int t{};
    std::cin >> t;
    if (t == 0) {
      int x{};
      int y{};
      int w{};
      int z{};
      std::cin >> x >> y >> w >> z;
      link_cut_tree.Cut(x, y);
      link_cut_tree.Link(w, z);
    } else if (t == 1) {
      int v{};
      int value{};
      std::cin >> v >> value;
      link_cut_tree.Apply(v, value);
    } else {
      int u{};
      int v{};
      std::cin >> u >> v;
      std::cout << link_cut_tree.GetPath(u, v) << "\n";
    }
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

