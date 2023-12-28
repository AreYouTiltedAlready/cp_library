#include <bits/stdc++.h>

namespace {

template <class Fun>
class y_combinator_result {
  Fun fun_;

 public:
  template <class T>
  explicit y_combinator_result(  // NOLINT(*-forwarding-reference-overload)
      T&& fun)
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

namespace splay_tree_traits {

// Inherit to implement subtree aggregate maintenance
template <typename S, auto Op>
class AggregateSplayTrait {
 protected:
  static_assert(
      std::is_convertible_v<decltype(Op), S (*)(S, S)> ||
          std::is_convertible_v<decltype(Op), S (*)(const S&, const S&)>,
      "Op must work as S(S, S)");

  template <typename SplayNode,
            typename std::enable_if_t<
                std::is_base_of_v<AggregateSplayTrait<S, Op>, SplayNode>,
                void>* = nullptr>
  void DoAggregatePull(const SplayNode* left, const SplayNode* right) {
    aggregate_value_ = node_value_;
    if (left != nullptr) {
      aggregate_value_ = Op(left->aggregate_value_, aggregate_value_);
    }
    if (right != nullptr) {
      aggregate_value_ = Op(aggregate_value_, right->aggregate_value_);
    }
  }
  explicit AggregateSplayTrait(const S& value)
      : node_value_(value), aggregate_value_(value) {}

  S node_value_;
  S aggregate_value_;
};

// Inherit to implement subtree size maintenance
class SizeSplayTrait {
 public:
  class Finder {
   public:
    template <
        typename SplayNode,
        typename std::enable_if_t<std::is_base_of_v<SizeSplayTrait, SplayNode>,
                                  void>* = nullptr>
    int operator()(SplayNode* node) const {
      int left_size = 0;
      if (auto* left = node->left(); left != nullptr) {
        left_size = left->size_;
      }
      if (left_size == k) {
        return 0;
      }
      if (left_size < k) {
        k -= left_size + 1;
        return 1;
      }
      return -1;
    }

    explicit Finder(int k) : k(k) {}

   private:
    mutable int k;
  };

  template <typename SplayNode,
            typename std::enable_if_t<
                std::is_base_of_v<SizeSplayTrait, SplayNode>, void>* = nullptr>
  friend SplayNode* FindKth(SplayNode* node, int k) {
    return Find(node, Finder(k));
  }

  template <typename SplayNode,
            typename std::enable_if_t<
                std::is_base_of_v<SizeSplayTrait, SplayNode>, void>* = nullptr>
  friend std::array<SplayNode*, 2> SplitKLeftmost(SplayNode* node, int k) {
    if (k <= 0) {
      return {{nullptr, node}};
    }
    if (k >= node->size_) {
      return {{node, nullptr}};
    }
    auto* aux = FindKth(node, k - 1);
    return SplitAfter(aux);
  }

 protected:
  SizeSplayTrait() : size_(1) {}

  template <typename SplayNode,
            typename std::enable_if_t<
                std::is_base_of_v<SizeSplayTrait, SplayNode>, void>* = nullptr>
  void DoSizePull(const SplayNode* left, const SplayNode* right) {
    size_ = 1;
    if (left != nullptr) {
      size_ += left->size_;
    }
    if (right != nullptr) {
      size_ += right->size_;
    }
  }

  int size_;
};

// Inherit to implement subtree reverse trait
class ReverseSplayTrait {
 public:
  template <
      typename SplayNode,
      typename std::enable_if_t<std::is_base_of_v<ReverseSplayTrait, SplayNode>,
                                void>* = nullptr>
  static void Reverse(SplayNode* node) {
    node->is_reversed_ ^= 1;
  }

 protected:
  template <
      typename SplayNode,
      typename std::enable_if_t<std::is_base_of_v<ReverseSplayTrait, SplayNode>,
                                void>* = nullptr>
  void DoReversePush(SplayNode** left, SplayNode** right) {
    if (is_reversed_ == 1) {
      std::swap(*left, *right);
      if (*left != nullptr) {
        Reverse(*left);
      }
      if (*right != nullptr) {
        Reverse(*right);
      }
      is_reversed_ = 0;
    }
  }

  ReverseSplayTrait() : is_reversed_(0) {}

  int is_reversed_;
};

}  // namespace splay_tree_traits

template <typename Policy>
class SplayTreeNode : public Policy {
 public:
  using internal_policy = Policy;
  using SplayPtr = SplayTreeNode<internal_policy>*;

  template <
      typename... Args,
      typename std::enable_if_t<
          std::is_constructible_v<internal_policy, Args...>, void>* = nullptr>
  explicit SplayTreeNode(Args&&... args)
      : internal_policy(std::forward<Args>(args)...) {
    left_ = nullptr;
    right_ = nullptr;
    parent_ = nullptr;
  }

  [[nodiscard]] SplayPtr left() const noexcept { return left_; }
  [[nodiscard]] SplayPtr right() const noexcept { return right_; }
  [[nodiscard]] SplayPtr parent() const noexcept { return parent_; }

  SplayPtr* mutable_left() noexcept { return &left_; }
  SplayPtr* mutable_right() noexcept { return &right_; }
  SplayPtr* mutable_parent() noexcept { return &parent_; }

  inline void Push() { internal_policy::PushPolicy(&left_, &right_); }
  inline void Pull() { internal_policy::PullPolicy(&left_, &right_); }

 private:
  void Rotate() {
    SplayPtr par = parent_;
    SplayPtr g = par->parent_;

    if (g != nullptr) {
      if (par == g->left_) {
        g->left_ = this;
      } else {
        g->right_ = this;
      }
    }

    par->Push();
    this->Push();
    if (this == par->left_) {
      par->left_ = this->right_;
      this->right_ = par;
      if (par->left_ != nullptr) {
        par->left_->parent_ = par;
      }
    } else {
      par->right_ = this->left_;
      this->left_ = par;
      if (par->right_ != nullptr) {
        par->right_->parent_ = par;
      }
    }

    this->parent_ = g;
    par->parent_ = this;
    par->Pull();
    this->Pull();
  }

  void Splay() {
    while (this->parent_ != nullptr) {
      if (SplayPtr g = parent_->parent_; g != nullptr) {
        if ((this == parent_->left_) == (parent_ == g->left_)) {
          parent_->Rotate();
        } else {
          this->Rotate();
        }
      }
      this->Rotate();
    }
    this->Push();
  }

  friend SplayPtr GetLeftmost(SplayPtr splay_ptr) {
    if (splay_ptr == nullptr) {
      return nullptr;
    }
    splay_ptr->Push();
    while (splay_ptr->left_ != nullptr) {
      splay_ptr = splay_ptr->left_;
      splay_ptr->Push();
    }
    splay_ptr->Splay();
    return splay_ptr;
  }

  friend SplayPtr GetRightmost(SplayPtr splay_ptr) {
    if (splay_ptr == nullptr) {
      return nullptr;
    }
    splay_ptr->Push();
    while (splay_ptr->right_ != nullptr) {
      splay_ptr = splay_ptr->right_;
      splay_ptr->Push();
    }
    splay_ptr->Splay();
    return splay_ptr;
  }

  friend std::array<SplayPtr, 2> SplitBefore(SplayPtr splay_ptr) {
    if (splay_ptr == nullptr) {
      return {{nullptr, nullptr}};
    }
    splay_ptr->Push();
    SplayPtr left = splay_ptr->left_;
    if (left != nullptr) {
      left->parent_ = nullptr;
      splay_ptr->left_ = nullptr;
      splay_ptr->Pull();
    }
    return {{left, splay_ptr}};
  }

  friend std::array<SplayPtr, 2> SplitAfter(SplayPtr splay_ptr) {
    if (splay_ptr == nullptr) {
      return {{nullptr, nullptr}};
    }
    splay_ptr->Push();
    SplayPtr right = splay_ptr->right_;
    if (right != nullptr) {
      right->parent_ = nullptr;
      splay_ptr->right_ = nullptr;
      splay_ptr->Pull();
    }
    return {{splay_ptr, right}};
  }

  template <typename Predicate>
  friend SplayPtr Find(SplayPtr splay_ptr, const Predicate& dir) {
    if (splay_ptr == nullptr) {
      return nullptr;
    }

    while (true) {
      splay_ptr->Push();
      if (int d = dir(splay_ptr); d != 0) {
        SplayPtr next = (d < 0 ? splay_ptr->left_ : splay_ptr->right_);
        if (next == nullptr) {
          break;
        }
        splay_ptr = next;
      } else {
        break;
      }
    }

    splay_ptr->Splay();
    return splay_ptr;
  }

  template <typename... Args>
  friend SplayPtr Merge(SplayPtr lhs, SplayPtr rhs, Args&&... args) {
    return Merge(Merge(lhs, rhs), args...);
  }

  friend SplayPtr Merge(SplayPtr lhs, SplayPtr rhs) {
    if (lhs == nullptr) {
      return rhs;
    }
    if (rhs == nullptr) {
      return lhs;
    }
    lhs->Push();
    rhs->Push();
    lhs = GetRightmost(lhs);
    lhs->right_ = rhs;
    rhs->parent_ = lhs;
    lhs->Pull();
    return lhs;
  }

  friend void Destroy(SplayPtr splay_ptr) {
    if (splay_ptr == nullptr) {
      return;
    }
    if (splay_ptr->left_ != nullptr) {
      splay_ptr->left_->parent_ = nullptr;
      Destroy(splay_ptr->left_);
      splay_ptr->left_ = nullptr;
    }
    if (splay_ptr->right_ != nullptr) {
      splay_ptr->right_->parent_ = nullptr;
      Destroy(splay_ptr->right_);
      splay_ptr->right_ = nullptr;
    }
    delete splay_ptr;
  }

  SplayPtr left_;
  SplayPtr right_;
  SplayPtr parent_;
};

template <typename SplayNode>
class SplayTree {
 public:
  using internal_policy = typename SplayNode::internal_policy;
  static constexpr bool kSupportsSizeTrait =
      std::is_base_of_v<splay_tree_traits::SizeSplayTrait, internal_policy>;

  SplayTree() : root_(nullptr) {}

  template <typename S,
            typename std::enable_if_t<
                std::is_constructible_v<internal_policy, S>, void>* = nullptr>
  explicit SplayTree(const std::vector<S>& values) : SplayTree() {
    if (values.empty()) {
      return;
    }
    for (const S& value : values) {
      root_ = Merge(root_, new SplayNode(value));
    }
  }

  SplayNode* CutRange(int left, int right) {
    auto [excess_left, desired_range, excess_right] = SearchRange(left, right);
    root_ = Merge(excess_left, excess_right);
    return desired_range;
  }

  void InsertRange(int pos, SplayNode* range) {
    auto [left, right] = SplitKLeftmost(root_, pos);
    root_ = Merge(left, range, right);
  }

  internal_policy GetRangePolicy(int left, int right) {
    auto [excess_left, desired_range, excess_right] = SearchRange(left, right);
    internal_policy result(*desired_range);  // NOLINT(*slicing)
    root_ = Merge(excess_left, desired_range, excess_right);
    return result;
  }

  void RemoveAt(int pos) {
    root_ = FindKth(root_, pos);
    SplayNode* left = *root_->mutable_left();
    SplayNode* right = *root_->mutable_right();
    if (left != nullptr) {
      *left->mutable_parent() = nullptr;
    }
    if (right != nullptr) {
      *right->mutable_parent() = nullptr;
    }
    delete root_;
    root_ = Merge(left, right);
  }

  template <
      typename... Args,
      typename std::enable_if_t<
          std::is_constructible_v<internal_policy, Args...>, void>* = nullptr>
  void InsertAt(int pos, Args&&... args) {
    auto [left, right] = SplitKLeftmost(root_, pos);
    root_ = Merge(left, new SplayNode(std::forward<Args>(args)...), right);
  }

  template <
      typename Callback,
      typename std::enable_if_t<
          std::is_invocable_r_v<void, Callback, SplayNode*>, void>* = nullptr>
  void ApplyOnRange(int left, int right, Callback&& cb) {
    auto [excess_left, desired_range, excess_right] = SearchRange(left, right);
    cb(desired_range);
    root_ = Merge(excess_left, desired_range, excess_right);
  }

  void Dump(std::ostream& os) { root_->Dump(os); }

  ~SplayTree() { Destroy(root_); }

 protected:
  std::array<SplayNode*, 3> SearchRange(int left, int right) {
    static_assert(kSupportsSizeTrait, "Policy must maintain subtree sizes");
    auto [excess_left, middle] = SplitKLeftmost(root_, left);
    auto [desired, excess_right] = SplitKLeftmost(middle, right - left);
    return {{excess_left, desired, excess_right}};
  }

  SplayNode* root_;
};

class Barrett {
 public:
  constexpr explicit Barrett(uint32_t mod)
      : mod_inverse(static_cast<uint64_t>(-1) / mod + 1), mod(mod) {}

  [[nodiscard]] uint32_t Product(uint32_t lhs, uint32_t rhs) const {
    return (*this)(static_cast<uint64_t>(lhs) * rhs);
  }

  [[nodiscard]] uint32_t operator()(uint64_t n) const {
    auto x = static_cast<uint64_t>(
        (static_cast<__uint128_t>(n) * mod_inverse) >> 64);
    uint64_t m = x * mod;
    return n - m + (n < m ? mod : 0);
  }

 private:
  uint64_t mod_inverse;
  uint32_t mod;
};

// Static modint class (for 32-bit compile-time modulos)
// Uses barrett reduction for multiplication. Quite fast in practice.
// Usage:
// using mint = StaticModint<998244353>; // whatever (1000000007, ...)
template <uint32_t kMod>
class StaticModint {
 public:
  static_assert(kMod < (1U << 30));

  StaticModint() : value_(0) {}

  template <typename T,
            typename std::enable_if_t<
                std::is_integral_v<T> && std::is_signed_v<T>, void>* = nullptr>
  StaticModint(T value)
      : value_(  // NOLINT(*-explicit-constructor)
            value % kMod) {
    if (value < 0) {
      value += kMod;
    }
  }

  template <typename T, typename std::enable_if_t<std::is_integral_v<T> &&
                                                      std::is_unsigned_v<T>,
                                                  void>* = nullptr>
  StaticModint(T value)
      : value_(  // NOLINT(*-explicit-constructor)
            value % kMod) {}

  using Mint = StaticModint<kMod>;

  static Mint Raw(unsigned value) {
    Mint result;
    result.value_ = value;
    return result;
  }

  [[nodiscard]] Mint operator+() const noexcept { return Mint(*this); }

  [[nodiscard]] Mint operator-() const noexcept { return Mint(kMod - value_); }

  [[nodiscard]] bool IsZero() const noexcept { return value_ == 0; }

  [[nodiscard]] explicit operator uint32_t() const { return value_; }

  [[nodiscard]] Mint Inverse() const noexcept {
    Mint m(*this);
    Mint result = Raw(1U);
#pragma GCC unroll(30)
    for (unsigned i = 0; i < 30; ++i) {
      if (((kMod - 2) >> i) % 2 == 1) {
        result *= m;
      }
      m *= m;
    }
    return result;
  }

  Mint& operator+=(const Mint& other) {
    if (value_ += other.value_; value_ >= kMod) {
      value_ -= kMod;
    }
    return *this;
  }

  Mint& operator-=(const Mint& other) {
    if (value_ += kMod - other.value_; value_ >= kMod) {
      value_ -= kMod;
    }
    return *this;
  }

  Mint& operator*=(const Mint& other) {
    value_ = barrett.Product(value_, other.value_);
    return *this;
  }

  Mint& operator/=(const Mint& other) {
    value_ = barrett.Product(value_, other.Inverse().value_);
    return *this;
  }

  friend Mint operator+(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) += rhs;
  }

  friend Mint operator-(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) -= rhs;
  }

  friend Mint operator*(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) *= rhs;
  }

  friend Mint operator/(const Mint& lhs, const Mint& rhs) {
    return Mint(lhs) /= rhs;
  }

  friend Mint Power(Mint m, uint64_t n) {
    Mint result = Raw(1U);
    while (n > 0) {
      if (n % 2 == 1) {
        result *= m;
      }
      m *= m;
      n /= 2;
    }
    return result;
  }

  friend bool operator==(const Mint& lhs, const Mint& rhs) {
    return lhs.value_ == rhs.value_;
  }

  friend std::istream& operator>>(std::istream& is, Mint& mint) {
    uint32_t value;
    is >> value;
    mint = Mint(value);
    return is;
  }

  friend std::ostream& operator<<(std::ostream& os, const Mint& mint) {
    return os << mint.value_;
  }

 private:
  static constexpr Barrett barrett = Barrett(kMod);
  uint32_t value_;
};

using Mint = StaticModint<998244353>;
inline Mint Op(const Mint& lhs, const Mint& rhs) { return lhs + rhs; }

class SplayTreePolicy : public splay_tree_traits::AggregateSplayTrait<Mint, Op>,
                        public splay_tree_traits::SizeSplayTrait,
                        public splay_tree_traits::ReverseSplayTrait {
 public:
  template <typename SplayNode,
            typename std::enable_if_t<
                std::is_base_of_v<SplayTreePolicy, SplayNode>, void>* = nullptr>
  friend void ApplyTag(SplayNode* node, Mint c, Mint d) {
    node->node_value_ *= c;
    node->node_value_ += d;
    node->aggregate_value_ *= c;
    node->aggregate_value_ += d * node->size_;
    node->a *= c;
    node->b *= c;
    node->b += d;
  }

  using SizeTrait = splay_tree_traits::SizeSplayTrait;
  using AggregateTrait = splay_tree_traits::AggregateSplayTrait<Mint, Op>;
  using ReverseTrait = splay_tree_traits::ReverseSplayTrait;

  explicit SplayTreePolicy(Mint value)
      : AggregateTrait(value), SizeTrait(), ReverseTrait(), a(1), b(0) {}

  [[nodiscard]] Mint value() const { return aggregate_value_; }

  template <typename SplayNode,
            typename std::enable_if_t<
                std::is_base_of_v<SplayTreePolicy, SplayNode>, void>* = nullptr>
  void PushPolicy(SplayNode** left, SplayNode** right) {
    DoReversePush(left, right);
    if (*left != nullptr) {
      ApplyTag(*left, a, b);
    }
    if (*right != nullptr) {
      ApplyTag(*right, a, b);
    }
    a = 1;
    b = 0;
  }

  template <typename SplayNode,
            typename std::enable_if_t<
                std::is_base_of_v<SplayTreePolicy, SplayNode>, void>* = nullptr>
  void PullPolicy(SplayNode** left, SplayNode** right) {
    DoSizePull(*left, *right);
    DoAggregatePull(*left, *right);
  }

  Mint a;
  Mint b;
};

// https://judge.yosupo.jp/problem/dynamic_sequence_range_affine_range_sum
// https://judge.yosupo.jp/submission/179998
void RunCase([[maybe_unused]] int testcase) {
  int n;
  int q;
  std::cin >> n >> q;

  std::vector<Mint> v(n);
  for (Mint& it : v) {
    std::cin >> it;
  }

  SplayTree<SplayTreeNode<SplayTreePolicy>> tree(v);
  while (q--) {
    int t;
    std::cin >> t;
    if (t == 0) {
      int id;
      Mint z;
      std::cin >> id >> z;
      tree.InsertAt(id, z);
    } else if (t == 1) {
      int id;
      std::cin >> id;
      tree.RemoveAt(id);
    } else {
      int left;
      int right;
      std::cin >> left >> right;
      if (t == 3) {
        Mint x;
        Mint y;
        std::cin >> x >> y;
        tree.ApplyOnRange(left, right,
                          [&](SplayTreeNode<SplayTreePolicy>* splay_node) {
                            ApplyTag(splay_node, x, y);
                          });
      } else if (t == 2) {
        tree.ApplyOnRange(left, right,
                          [](SplayTreeNode<SplayTreePolicy>* splay_node) {
                            SplayTreePolicy::ReverseTrait::Reverse(splay_node);
                          });
      } else {
        std::cout << tree.GetRangePolicy(left, right).value() << "\n";
      }
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
