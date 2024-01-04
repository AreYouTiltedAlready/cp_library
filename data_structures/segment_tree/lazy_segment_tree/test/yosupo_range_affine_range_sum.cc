// Problem: https://judge.yosupo.jp/problem/range_affine_range_sum
// Submission: https://judge.yosupo.jp/submission/181563

#include <algorithm>
#include <functional>
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

namespace ds {
namespace lazy_segment_tree {

namespace internal {

template <typename T>
constexpr bool has_binary_plus = requires(const T& lhs, const T& rhs) {
  { lhs + rhs } -> std::same_as<T>;
};

template <typename T>
concept monoid = std::is_default_constructible_v<T> && has_binary_plus<T>;

template <typename T, typename U>
concept lazy_tag = has_binary_plus<T> && requires(U& u, const T& t, int z) {
  { u.Apply(t, z) } -> std::same_as<void>;
};

}  // namespace internal

template <internal::monoid S, internal::lazy_tag<S> T>
class LazySegmentTree {
 public:
  explicit LazySegmentTree(int n)
      : tree_(std::bit_ceil(static_cast<uint32_t>(n)) * 2),
        tag_(std::bit_ceil(static_cast<uint32_t>(n))),
        n_(n),
        size_(static_cast<int>(std::bit_ceil(static_cast<uint32_t>(n)))),
        log_(std::bit_width(static_cast<uint32_t>(size_) - 1)) {}

  template <typename U>
    requires std::is_assignable_v<S&, U>
  explicit LazySegmentTree(int n, const U& value) : LazySegmentTree(n) {
    std::ranges::fill(tree_, value);
    for (int i = size_ - 1; i > 0; --i) {
      Pull(i);
    }
  }

  template <typename U>
    requires std::is_assignable_v<S&, U>
  explicit LazySegmentTree(const std::vector<U>& values)
      : LazySegmentTree(static_cast<int>(values.size())) {
    std::ranges::copy(values, tree_.begin() + size_);
    for (int i = size_ - 1; i > 0; --i) {
      Pull(i);
    }
  }

  template <typename U>
    requires std::is_assignable_v<S&, U>
  void Set(int pos, const U& value) {
    pos += size_;
    for (int i = log_; i > 0; --i) {
      Push(pos >> i);
    }
    tree_[pos] = value;
    for (int i = 1; i <= log_; ++i) {
      Pull(pos >> i);
    }
  }

  void Apply(int pos, const T& tag) {
    pos += size_;
    for (int i = log_; i > 0; --i) {
      Push(pos >> i);
    }
    tree_[pos].Apply(tag, 1);
    for (int i = 1; i <= log_; ++i) {
      Pull(pos >> i);
    }
  }

  [[nodiscard]] S Get() const { return tree_[1]; }

  S Get(int pos) {
    pos += size_;
    for (int i = log_; i > 0; --i) {
      Push(pos >> i);
    }
    return tree_[pos];
  }

  void Apply(int first, int last, const T& tag) {
    first += size_;
    last += size_;
    DoPush(first, last);
    {
      int first_copy = first;
      int last_copy = last;
      while (first < last) {
        if ((first & 1) == 1) {
          ApplyAll(first++, tag);
        }
        if ((last & 1) == 1) {
          ApplyAll(--last, tag);
        }
        first >>= 1;
        last >>= 1;
      }
      first = first_copy;
      last = last_copy;
    }
    DoPull(first, last);
  }

  S Get(int first, int last) {
    first += size_;
    last += size_;
    DoPush(first, last);
    S res_left{};
    S res_right{};
    while (first < last) {
      if ((first & 1) == 1) {
        res_left = res_left + tree_[first++];
      }
      if ((last & 1) == 1) {
        res_right = tree_[--last] + res_right;
      }
      first >>= 1;
      last >>= 1;
    }
    return res_left + res_right;
  }

  using Predicate = std::function<bool(const S&)>;

  [[nodiscard]] int FindFirst(int first, int last, const Predicate& pred) {
    first += size_;
    last += size_;
    DoPush(first, last);
    int first_copy = first;
    int last_copy = last;
    while (first_copy < last_copy) {
      if ((first_copy & 1) == 1) {
        if (pred(tree_[first_copy])) {
          return Descent(first_copy, pred, DescentDirection::kLeft);
        }
        first_copy += 1;
      }

      first_copy >>= 1;
      last_copy >>= 1;
    }

    int height = std::bit_width(static_cast<uint32_t>(last_copy));
    for (int i = log_ - height; i >= 0; --i) {
      last_copy <<= 1;
      if (((last >> i) & 1) == 1) {
        if (pred(tree_[last_copy])) {
          return Descent(last_copy, pred, DescentDirection::kLeft);
        }
        last_copy += 1;
      }
    }

    return -1;
  }

  [[nodiscard]] int FindLast(int first, int last, const Predicate& pred) {
    first += size_;
    last += size_;
    DoPush(first, last);

    uint32_t mask = 1;
    int first_copy = first;
    int last_copy = last;
    while (first_copy < last_copy) {
      mask <<= 1;
      if ((first_copy & 1) == 1) {
        mask += 1;
        first_copy += 1;
      }

      if ((last_copy & 1) == 1) {
        last_copy -= 1;
        if (pred(tree_[last_copy])) {
          return Descent(last_copy, pred, DescentDirection::kRight);
        }
      }

      first_copy >>= 1;
      last_copy >>= 1;
    }

    while (!std::has_single_bit(mask)) {
      int z = std::countr_zero(mask) + 1;
      mask >>= z;
      first_copy <<= z;
      first_copy -= 1;
      if (pred(tree_[first_copy])) {
        return Descent(first_copy, pred, DescentDirection::kRight);
      }
    }

    return -1;
  }

 private:
  enum class DescentDirection {
    kLeft = 0,
    kRight = 1,
  };

  [[nodiscard]] int Descent(int k, const Predicate& pred,
                            DescentDirection direction) {
    while (k < size_) {
      Push(k);
      k <<= 1;
      k ^= static_cast<int>(direction);
      k ^= !pred(tree_[k]);
    }
    return k - size_;
  }

  void ApplyAll(int k, const T& tag) {
    int length = size_ >> (std::bit_width(static_cast<uint32_t>(k)) - 1);
    tree_[k].Apply(tag, length);
    if (length != 1) [[likely]] {
      tag_[k] = tag_[k] + tag;
    }
  }

  void DoPull(int first, int last) {
    const int last_common_bit = std::max(
        1, static_cast<int>(
               std::bit_width(static_cast<uint32_t>(first ^ (last - 1)))));
    int first_z = std::max(1, std::countr_zero(static_cast<uint32_t>(first)));
    for (int i = first_z; i < last_common_bit; ++i) {
      Pull(first >> i);
    }
    int last_z = std::max(1, std::countr_zero(static_cast<uint32_t>(last)));
    for (int i = last_z; i < last_common_bit; ++i) {
      Pull((last - 1) >> i);
    }
    for (int i = last_common_bit; i <= log_; ++i) {
      Pull(first >> i);
    }
  }

  void DoPush(int first, int last) {
    const int last_common_bit = std::max(
        1, static_cast<int>(
               std::bit_width(static_cast<uint32_t>(first ^ (last - 1)))));
    for (int i = log_; i >= last_common_bit; --i) {
      Push(first >> i);
    }
    int first_z = std::max(1, std::countr_zero(static_cast<uint32_t>(first)));
    for (int i = last_common_bit - 1; i >= first_z; --i) {
      Push(first >> i);
    }
    int last_z = std::max(1, std::countr_zero(static_cast<uint32_t>(last)));
    for (int i = last_common_bit - 1; i >= last_z; --i) {
      Push((last - 1) >> i);
    }
  }

  void Push(int k) {
    ApplyAll(k << 1, tag_[k]);
    ApplyAll(k << 1 | 1, tag_[k]);
    tag_[k] = T();
  }

  inline void Pull(int k) {
    tree_[k] = tree_[k << 1] + tree_[k << 1 | 1];
    int length = size_ >> (std::bit_width(static_cast<uint32_t>(k)) - 1);
    tree_[k].Apply(tag_[k], length);
  }

  std::vector<S> tree_;
  std::vector<T> tag_;
  const int n_;
  const int size_;
  const int log_;
};

}  // namespace lazy_segment_tree

}  // namespace ds

namespace modular {

namespace internal {
namespace integral_type_traits {

template <typename T>
struct integral_promotion {
  using type = T;
};

template <>
struct integral_promotion<int> {
  using type = int64_t;
};

template <>
struct integral_promotion<uint32_t> {
  using type = uint64_t;
};

template <>
struct integral_promotion<int64_t> {
  using type = __int128_t;
};

template <>
struct integral_promotion<uint64_t> {
  using type = __uint128_t;
};

}  // namespace integral_type_traits

template <typename T>
using integral_promotion_t =
    typename integral_type_traits::integral_promotion<T>::type;

namespace montgomery {

template <typename T>
class MontgomerySpace {
 public:
  using signed_t = std::make_signed_t<T>;
  using unsigned_t = std::make_unsigned_t<T>;
  using promoted_t = integral_promotion_t<unsigned_t>;
  static constexpr int kTBitWidth = std::numeric_limits<unsigned_t>::digits;

  constexpr explicit MontgomerySpace(T mod)
      : mod_(mod),
        r_square_((static_cast<promoted_t>(1) << kTBitWidth) % mod),
        mod_inverse_(1) {
    r_square_ = static_cast<promoted_t>(r_square_) * r_square_ % mod;
    for (int i = 0; i < 6; ++i) {
      mod_inverse_ *= static_cast<unsigned_t>(2) - mod_ * mod_inverse_;
    }
  }

  [[nodiscard]] constexpr unsigned_t Sum(unsigned_t lhs, unsigned_t rhs) const {
    if (static_cast<signed_t>(lhs += rhs - (mod_ << 1)) < 0) {
      lhs += mod_ << 1;
    }
    return lhs;
  }

  [[nodiscard]] constexpr unsigned_t Difference(unsigned_t lhs,
                                                unsigned_t rhs) const {
    if (static_cast<signed_t>(lhs -= rhs) < 0) {
      lhs += mod_ << 1;
    }
    return lhs;
  }

  [[nodiscard]] constexpr unsigned_t Product(unsigned_t lhs,
                                             unsigned_t rhs) const {
    return Reduce(static_cast<promoted_t>(lhs) * rhs);
  }

  [[nodiscard]] constexpr bool AreEqual(unsigned_t lhs, unsigned_t rhs) const {
    return (lhs < mod_ ? lhs : lhs - mod_) == (rhs < mod_ ? rhs : rhs - mod_);
  }

  [[nodiscard]] constexpr unsigned_t Power(unsigned_t x, uint64_t n) const {
    unsigned_t result = Transform(1U);
    while (n > 0) {
      if (n % 2 == 1) {
        result = Product(result, x);
      }
      x = Product(x, x);
      n /= 2;
    }
    return result;
  }

  // returns a value congruent to mod in range [0, 2 * mod)
  [[nodiscard]] constexpr unsigned_t Reduce(promoted_t n) const {
    unsigned_t q = static_cast<unsigned_t>(n) * mod_inverse_;
    unsigned_t m = (static_cast<promoted_t>(q) * mod_) >> kTBitWidth;
    return (n >> kTBitWidth) + mod_ - m;
  }

  [[nodiscard]] constexpr unsigned_t Transform(unsigned_t n) const {
    return Product(n, r_square_);
  }

 private:
  unsigned_t mod_;
  unsigned_t r_square_;
  unsigned_t mod_inverse_;
};

}  // namespace montgomery

struct ModIntBase {};

template <typename T>
using is_modint = std::is_base_of<ModIntBase, T>;

template <typename T>
using is_modint_t = typename std::is_base_of<ModIntBase, T>::type;

template <typename T, T kMod>
class LazyMontgomeryModInt : public internal::ModIntBase {
  using signed_t = std::make_signed_t<T>;
  using unsigned_t = std::make_unsigned_t<T>;

  static_assert(kMod > 0, "Mod must be positive");
  static_assert(kMod % 2 == 1, "Mod must be odd");
  static_assert((kMod < (static_cast<uint64_t>(1)
                         << (std::numeric_limits<signed_t>::digits - 1))),
                "Mod is too large");

 public:
  using Mint = LazyMontgomeryModInt<T, kMod>;

  constexpr LazyMontgomeryModInt() : value_(0) {}

  template <typename U,
            typename std::enable_if_t<std::is_signed_v<U>>* = nullptr>
  constexpr LazyMontgomeryModInt(U n)  // NOLINT(*explicit-constructor*)
      : value_(0) {
    if (n %= kMod; n < 0) {
      n += kMod;
    }
    value_ = space_.Transform(static_cast<unsigned_t>(n));
  }

  template <typename U,
            typename std::enable_if_t<std::is_unsigned_v<U>>* = nullptr>
  constexpr LazyMontgomeryModInt(U n)  // NOLINT(*explicit-constructor*)
      : value_(space_.Transform(static_cast<unsigned_t>(n % kMod))) {}

  static constexpr unsigned_t UMod() { return kUMod_; }

  static constexpr Mint Raw(unsigned_t value) {
    Mint result;
    result.value_ = space_.Transform(value);
    return result;
  }

  [[nodiscard]] constexpr Mint operator+() const noexcept { return *this; }

  [[nodiscard]] constexpr Mint operator-() const noexcept {
    return Mint() - *this;
  }

  template <typename U,
            typename std::enable_if_t<std::is_integral_v<U>>* = nullptr>
  [[nodiscard]] constexpr explicit operator U() const {
    unsigned_t q = space_.Reduce(value_);
    return static_cast<T>(q < kUMod_ ? q : q - kUMod_);
  }

  [[nodiscard]] constexpr Mint Inverse() const noexcept {
    return Power(*this, kUMod_ - 2);
  }

  constexpr Mint& operator+=(const Mint& other) noexcept {
    value_ = space_.Sum(value_, other.value_);
    return *this;
  }

  constexpr Mint& operator-=(const Mint& other) noexcept {
    value_ = space_.Difference(value_, other.value_);
    return *this;
  }

  constexpr Mint& operator*=(const Mint& other) noexcept {
    value_ = space_.Product(value_, other.value_);
    return *this;
  }

  constexpr Mint& operator/=(const Mint& other) noexcept {
    value_ = space_.Product(value_, other.Inverse().value_);
    return *this;
  }

  constexpr friend Mint operator+(const Mint& lhs, const Mint& rhs) {
    Mint res{};
    res.value_ = space_.Sum(lhs.value_, rhs.value_);
    return res;
  }

  constexpr friend Mint operator-(const Mint& lhs, const Mint& rhs) {
    Mint res{};
    res.value_ = space_.Difference(lhs.value_, rhs.value_);
    return res;
  }

  constexpr friend Mint operator*(const Mint& lhs, const Mint& rhs) {
    Mint res{};
    res.value_ = space_.Product(lhs.value_, rhs.value_);
    return res;
  }

  constexpr friend Mint operator/(const Mint& lhs, const Mint& rhs) {
    Mint res{};
    res.value_ = space_.Product(lhs.value_, rhs.Inverse().value_);
    return res;
  }

  constexpr friend bool operator==(const Mint& lhs, const Mint& rhs) {
    return space_.AreEqual(lhs.value_, rhs.value_);
  }

  constexpr friend bool operator!=(const Mint& lhs, const Mint& rhs) {
    return !space_.AreEqual(lhs.value_, rhs.value_);
  }

  constexpr friend Mint Power(Mint mint, uint64_t n) noexcept {
    Mint res;
    res.value_ = space_.Power(mint.value_, n);
    return res;
  }

  friend std::istream& operator>>(std::istream& istream, Mint& mint) {
    unsigned_t value;
    istream >> value;
    mint = Mint::Raw(value);
    return istream;
  }

  friend std::ostream& operator<<(std::ostream& ostream, const Mint& mint) {
    return ostream << static_cast<signed_t>(mint);
  }

 private:
  static constexpr unsigned_t kUMod_ = kMod;
  static constexpr internal::montgomery::MontgomerySpace<unsigned_t> space_ =
      internal::montgomery::MontgomerySpace<unsigned_t>(kUMod_);

  unsigned_t value_;
};

}  // namespace internal

template <int kMod>
using MInt = internal::LazyMontgomeryModInt<int, kMod>;

template <int64_t kMod>
using MLong = internal::LazyMontgomeryModInt<int64_t, kMod>;

}  // namespace modular

using Mint = modular::MInt<998244353>;

struct Tag {
  Tag() : k(Mint::Raw(1)), b() {}
  constexpr Tag(Mint k, Mint b) : k(k), b(b) {}

  friend Tag operator+(const Tag& lhs, const Tag& rhs) {
    return {lhs.k * rhs.k, rhs.k * lhs.b + rhs.b};
  }

  Mint k;
  Mint b;
};

struct S {
  S() : value() {}
  S(Mint value) : value(value) {}

  friend S operator+(const S& lhs, const S& rhs) {
    return lhs.value + rhs.value;
  }

  void Apply(const Tag& tag, int length) {
    value *= tag.k;
    value += tag.b * length;
  }

  Mint value;
};

void RunCase([[maybe_unused]] int testcase) {
  int n;
  int q;
  std::cin >> n >> q;

  std::vector<Mint> v(n);
  for (Mint& it : v) {
    std::cin >> it;
  }

  ds::lazy_segment_tree::LazySegmentTree<S, Tag> segment_tree(v);
  while (q--) {
    int t;
    int first;
    int last;
    std::cin >> t >> first >> last;
    if (t == 0) {
      Mint k;
      Mint b;
      std::cin >> k >> b;
      segment_tree.Apply(first, last, Tag(k, b));
    } else {
      std::cout << segment_tree.Get(first, last).value << "\n";
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
