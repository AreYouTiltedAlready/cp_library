#include <immintrin.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <utility>

namespace bitset {
namespace internal {

constexpr __m256i kZeroes256 = {0, 0, 0, 0};
constexpr __m256i kOnes256 = {
    static_cast<int64_t>(-1), static_cast<int64_t>(-1),
    static_cast<int64_t>(-1), static_cast<int64_t>(-1)};

constexpr inline uint64_t Mask(int bits) {
  return bits == 64 ? std::numeric_limits<uint64_t>::max()
                    : (static_cast<uint64_t>(1) << bits) - 1;
}

constexpr inline uint64_t Bit(int bit) {
  return static_cast<uint64_t>(1) << bit;
}

constexpr inline int Ceil256(int n) {
  return static_cast<int>((n + Mask(8))) >> 8;
}

constexpr inline int Blocks(int n) { return Ceil256(n) << 2; }

uint64_t* Allocate(int n) {
  int blocks = Blocks(n);
  auto* res = std::assume_aligned<32>(
      new (static_cast<std::align_val_t>(32))  // NOLINT(*owning-memory*)
      uint64_t[blocks]);
  for (int i = 0; i < blocks; i += 4) {
    _mm256_stream_si256(reinterpret_cast<__m256i*>(res + i), kZeroes256);
  }
  return res;
}

void Deallocate(uint64_t* data) {
  ::operator delete[](data, static_cast<std::align_val_t>(32));
}

}  // namespace internal

class avx_bitset {
 public:
  avx_bitset() noexcept
      : data_(nullptr), n_(0), blocks_256_(0), blocks_64_(0) {}

  explicit avx_bitset(int n)
      : data_(internal::Allocate(n)),
        n_(n),
        blocks_256_(n >> 8),
        blocks_64_(n >> 6) {}

  static constexpr uint64_t kZeroMask = 0;
  static constexpr uint64_t kOneMask = internal::Mask(64);

  avx_bitset(const avx_bitset& other) : avx_bitset() { CopyFrom(other); }

  avx_bitset(avx_bitset&& other) noexcept : avx_bitset() {
    MoveFrom(std::move(other));
  }

  avx_bitset& operator=(const avx_bitset& other) {
    if (this == &other) {
      return *this;
    }
    internal::Deallocate(data_);
    CopyFrom(other);
    return *this;
  }

  avx_bitset& operator=(avx_bitset&& other) noexcept {
    if (this == &other) {
      return *this;
    }
    internal::Deallocate(data_);
    MoveFrom(std::move(other));
    return *this;
  }

  [[nodiscard]] bool test(int pos) const {
    return (data_[pos >> 6] & internal::Bit(pos & 63)) != 0;
  }

  void set(int pos) { data_[pos >> 6] |= internal::Bit(pos & 63); }

  void flip(int pos) { data_[pos >> 6] ^= internal::Bit(pos & 63); }

  void reset(int pos) { data_[pos >> 6] &= kOneMask ^ internal::Bit(pos & 63); }

  [[nodiscard]] int count() const {
    int result = 0;
    for (int i = 0; i < blocks_64_; ++i) {
      result += std::popcount(*(data_ + i));
    }
    result +=
        ((n_ & 63) == 0
             ? 0
             : std::popcount(data_[blocks_64_] & internal::Mask(n_ & 63)));
    return result;
  }

  [[nodiscard]] bool none() const { return !any(); }

  [[nodiscard]] bool any() const {
    for (int i = 0; i < blocks_256_; ++i) {
      __m256i block = _mm256_stream_load_si256(reinterpret_cast<const __m256i*>(
          data_ + (static_cast<ptrdiff_t>(i) * 4)));
      if (_mm256_testz_si256(block, block) == 0) {
        return true;
      }
    }
    for (int i = blocks_256_ << 2; i < blocks_64_; ++i) {
      if (data_[i] != kZeroMask) {
        return true;
      }
    }
    return (n_ & 63) != 0 && (data_[blocks_64_] & internal::Mask(n_ & 63)) != 0;
  }

  [[nodiscard]] bool all() const {
    for (int i = 0; i < blocks_256_; ++i) {
      __m256i flipped_block = _mm256_xor_si256(
          _mm256_stream_load_si256(reinterpret_cast<const __m256i*>(
              data_ + (static_cast<ptrdiff_t>(i) * 4))),
          internal::kOnes256);
      if (_mm256_testz_si256(flipped_block, flipped_block) == 0) {
        return false;
      }
    }
    for (int i = blocks_256_ << 2; i < blocks_64_; ++i) {
      if (data_[i] != kOneMask) {
        return false;
      }
    }
    return (n_ & 63) == 0 || (data_[blocks_64_] & internal::Mask(n_ & 63)) ==
                                 internal::Mask(n_ & 63);
  }

  // Requires (*this).n_ == other.n_
  avx_bitset& operator|=(const avx_bitset& other) {
    for (int i = 0; i < blocks_256_; ++i) {
      __m256i block =
          _mm256_or_si256(_mm256_load_si256(reinterpret_cast<const __m256i*>(
                              data_ + (static_cast<ptrdiff_t>(i) * 4))),
                          _mm256_load_si256(reinterpret_cast<const __m256i*>(
                              other.data_ + (static_cast<ptrdiff_t>(i) * 4))));
      _mm256_stream_si256(
          reinterpret_cast<__m256i*>(data_ + (static_cast<ptrdiff_t>(i) * 4)),
          block);
    }
    for (int i = blocks_256_ << 2; i < blocks_64_; ++i) {
      data_[i] |= other.data_[i];
    }
    if ((n_ & 63) != 0) {
      data_[blocks_64_] |= other.data_[blocks_64_];
    }
    return *this;
  }

  // Requires (*this).n_ == other.n_
  avx_bitset& operator&=(const avx_bitset& other) {
    for (int i = 0; i < blocks_256_; ++i) {
      __m256i block =
          _mm256_and_si256(_mm256_load_si256(reinterpret_cast<const __m256i*>(
                               data_ + (static_cast<ptrdiff_t>(i) * 4))),
                           _mm256_load_si256(reinterpret_cast<const __m256i*>(
                               other.data_ + (static_cast<ptrdiff_t>(i) * 4))));
      _mm256_stream_si256(
          reinterpret_cast<__m256i*>(data_ + (static_cast<ptrdiff_t>(i) * 4)),
          block);
    }
    for (int i = blocks_256_ << 2; i < blocks_64_; ++i) {
      data_[i] &= other.data_[i];
    }
    if ((n_ & 63) != 0) {
      data_[blocks_64_] &= other.data_[blocks_64_];
    }
    return *this;
  }

  // Requires (*this).n_ == other.n_
  avx_bitset& operator^=(const avx_bitset& other) {
    for (int i = 0; i < blocks_256_; ++i) {
      __m256i block =
          _mm256_xor_si256(_mm256_load_si256(reinterpret_cast<const __m256i*>(
                               data_ + (static_cast<ptrdiff_t>(i) * 4))),
                           _mm256_load_si256(reinterpret_cast<const __m256i*>(
                               other.data_ + (static_cast<ptrdiff_t>(i) * 4))));
      _mm256_stream_si256(
          reinterpret_cast<__m256i*>(data_ + (static_cast<ptrdiff_t>(i) * 4)),
          block);
    }
    for (int i = blocks_256_ << 2; i < blocks_64_; ++i) {
      data_[i] ^= other.data_[i];
    }
    if ((n_ & 63) != 0) {
      data_[blocks_64_] ^= other.data_[blocks_64_];
    }
    return *this;
  }

  friend bool operator==(const avx_bitset& lhs, const avx_bitset& rhs) {
    // This check is needed because we read data non-temporal, thus
    // caching it in case of equal addresses would be better
    if (&lhs == &rhs) {
      return true;
    }

    if (lhs.n_ != rhs.n_) {
      return false;
    }

    for (int i = 0; i < lhs.blocks_256_; ++i) {
      __m256i xor_block = _mm256_xor_si256(
          _mm256_stream_load_si256(reinterpret_cast<const __m256i*>(
              lhs.data_ + (static_cast<ptrdiff_t>(i) * 4))),
          _mm256_stream_load_si256(reinterpret_cast<const __m256i*>(
              rhs.data_ + (static_cast<ptrdiff_t>(i) * 4))));
      if (_mm256_testz_si256(xor_block, xor_block) == 0) {
        return false;
      }
    }

    for (int i = 0; i < lhs.blocks_64_; ++i) {
      if (lhs.data_[i] != rhs.data_[i]) {
        return false;
      }
    }

    return (lhs.n_ & 63) == 0 ||
           lhs.data_[lhs.blocks_64_] == rhs.data_[rhs.blocks_64_];
  }

  friend bool operator!=(const avx_bitset& lhs, const avx_bitset& rhs) {
    return !(lhs == rhs);
  }

  friend avx_bitset operator|(const avx_bitset& lhs, const avx_bitset& rhs) {
    avx_bitset res(lhs);
    res |= rhs;
    return res;
  }

  friend avx_bitset operator&(const avx_bitset& lhs, const avx_bitset& rhs) {
    avx_bitset res(lhs);
    res &= rhs;
    return res;
  }

  friend avx_bitset operator^(const avx_bitset& lhs, const avx_bitset& rhs) {
    avx_bitset res(lhs);
    res ^= rhs;
    return res;
  }

  avx_bitset& operator<<=(int shift) {
    if (shift == 0) {
      return *this;
    }
    if (shift >= n_) {
      reset();
      return *this;
    }
    if (int shift_blocks = shift >> 6; shift_blocks != 0) {
      auto* new_data = internal::Allocate(n_);
      std::memcpy(new_data + shift_blocks, data_,
                  (internal::Blocks(n_) - shift_blocks) * sizeof(uint64_t));
      std::swap(data_, new_data);
      internal::Deallocate(data_);
    }
    if ((shift & 63) != 0) {
      ShiftRightSmall(shift & 63);
    }
    RemoveTrash();
    return *this;
  }

  avx_bitset operator<<(int shift) const {
    if (shift == 0) {
      return avx_bitset(*this);
    }
    if (shift >= n_) {
      return avx_bitset(n_);
    }
    auto* new_data = internal::Allocate(n_);
    int shift_blocks = shift >> 6;
    std::memcpy(new_data + shift_blocks, data_,
                (internal::Blocks(n_) - shift_blocks) * sizeof(uint64_t));
    avx_bitset result(new_data, n_);
    if ((shift & 63) != 0) {
      result.ShiftRightSmall(shift & 63);
    }
    result.RemoveTrash();
    return result;
  }

  avx_bitset& operator>>=(int shift) {
    if (shift == 0) {
      return *this;
    }
    if (shift >= n_) {
      reset();
      return *this;
    }
    if (int shift_blocks = shift >> 6; shift_blocks != 0) {
      auto* new_data = internal::Allocate(n_);
      std::memcpy(new_data, data_ + shift_blocks,
                  (internal::Blocks(n_) - shift_blocks) * sizeof(uint64_t));
      std::swap(data_, new_data);
      internal::Deallocate(data_);
    }
    if ((shift & 63) != 0) {
      ShiftLeftSmall(shift & 63);
    }
    return *this;
  }

  avx_bitset operator>>(int shift) const {
    if (shift == 0) {
      return avx_bitset(*this);
    }
    if (shift >= n_) {
      return avx_bitset(n_);
    }
    auto* new_data = internal::Allocate(n_);
    int shift_blocks = shift >> 6;
    std::memcpy(new_data, data_ + shift_blocks,
                (internal::Blocks(n_) - shift_blocks) * sizeof(uint64_t));
    avx_bitset result(new_data, n_);
    if ((shift & 63) != 0) {
      result.ShiftLeftSmall(shift & 63);
    }
    return result;
  }

  [[nodiscard]] int find_first_one(int first, int last) const {
    if ((first >> 6) == (last >> 6)) {
      while (first < last) {
        if (test(first)) {
          return first;
        }
        first += 1;
      }
      return -1;
    }

    if (int r = first & 63; r != 0) {
      if (uint64_t suffix = data_[first >> 6] >> r; suffix != 0) {
        return first + std::countr_zero(suffix);
      }
      first += 64 - r;
    }

    int first_block = first >> 6;
    int last_block = last >> 6;
    while (first_block < last_block && (first_block & 3) != 0) {
      if (data_[first_block] != 0) {
        return (first_block << 6) + std::countr_zero(data_[first_block]);
      }
      first_block += 1;
    }

    while (first_block + 4 <= last_block) {
      __m256i block = _mm256_stream_load_si256(
          reinterpret_cast<const __m256i*>(data_ + first_block));
      if (_mm256_testz_si256(block, block) == 0) {
        while (data_[first_block] == 0) {
          first_block += 1;
        }
        return (first_block << 6) + std::countr_zero(data_[first_block]);
      }
      first_block += 4;
    }

    while (first_block < last_block) {
      if (data_[first_block] != 0) {
        return (first_block << 6) + std::countr_zero(data_[first_block]);
      }
      first_block += 1;
    }

    if (int r = last & 63; r != 0) {
      if (uint64_t prefix = data_[last >> 6] & internal::Mask(r); prefix != 0) {
        return last - r + std::countr_zero(prefix);
      }
    }

    return -1;
  }

  [[nodiscard]] int find_last_one(int first, int last) const {
    if ((first >> 6) == (last >> 6)) {
      while (last > first) {
        if (test(--last)) {
          return last;
        }
      }
      return -1;
    }

    int first_block = first >> 6;
    int last_block = last >> 6;
    if (int r = last & 63; r != 0) {
      if (uint64_t prefix = data_[last_block] & internal::Mask(r);
          prefix != 0) {
        return (last_block << 6) - std::countl_zero(prefix) + 63;
      }
      last -= r;
    }

    while (first_block < last_block && (last_block & 3) != 0) {
      if (data_[--last_block] != kZeroMask) {
        return (last_block << 6) - std::countl_zero(data_[last_block]) + 63;
      }
    }

    while (last_block - 4 > first_block) {
      __m256i block = _mm256_stream_load_si256(
          reinterpret_cast<const __m256i*>(data_ + last_block - 4));
      if (_mm256_testz_si256(block, block) == 0) {
        while (data_[last_block - 1] == kZeroMask) {
          last_block -= 1;
        }
        last_block -= 1;
        return (last_block << 6) - std::countl_zero(data_[last_block]) + 63;
      }
      last_block -= 4;
    }

    while (last_block - 1 > first_block) {
      if (data_[--last_block] != kZeroMask) {
        return (last_block << 6) - std::countl_zero(data_[last_block]) + 63;
      }
    }

    int r = first & 63;
    if (uint64_t suffix = data_[first_block] >> r; suffix != 0) {
      return (first_block << 6) - std::countl_zero(suffix) + 63;
    }

    return -1;
  }

  avx_bitset& flip() {
    for (int i = 0; i < blocks_256_; ++i) {
      __m256i flipped_block = _mm256_xor_si256(
          _mm256_load_si256(reinterpret_cast<const __m256i*>(data_ + (i << 2))),
          internal::kOnes256);
      _mm256_stream_si256(reinterpret_cast<__m256i*>(data_ + (i << 2)),
                          flipped_block);
    }
    for (int i = blocks_256_ << 2; i < blocks_64_; ++i) {
      data_[i] ^= kOneMask;
    }
    if (int r = n_ & 63; r != 0) {
      data_[blocks_64_] ^= internal::Mask(r);
    }
    return *this;
  }

  avx_bitset& reset() {
    for (int i = 0; i < blocks_256_; ++i) {
      _mm256_stream_si256(reinterpret_cast<__m256i*>(data_ + (i << 2)),
                          internal::kZeroes256);
    }
    for (int i = blocks_256_ << 2; i < blocks_64_; ++i) {
      data_[i] = kZeroMask;
    }
    if (int r = n_ & 63; r != 0) {
      data_[blocks_64_] = kZeroMask;
    }
    return *this;
  }

  avx_bitset& set() {
    for (int i = 0; i < blocks_256_; ++i) {
      _mm256_stream_si256(reinterpret_cast<__m256i*>(data_ + (i << 2)),
                          internal::kOnes256);
    }
    for (int i = blocks_256_ << 2; i < blocks_64_; ++i) {
      data_[i] = kOneMask;
    }
    if (int r = n_ & 63; r != 0) {
      data_[blocks_64_] |= internal::Mask(r);
    }
    return *this;
  }

  explicit operator std::string() const {
    std::string result(n_, '0');
    for (int i = 0; i < blocks_64_; ++i) {
      for (int j = 0; j < 64; ++j) {
        if ((data_[i] & internal::Bit(j)) != 0) {
          result[(i << 6) + j] = '1';
        }
      }
    }
    for (int i = 0; i < (n_ & 63); ++i) {
      if ((data_[blocks_64_] & internal::Bit(i)) != 0) {
        result[(blocks_64_ << 6) + i] = '1';
      }
    }
    std::reverse(result.begin(), result.end());
    return result;
  }

  [[nodiscard]] std::string to_string() const {
    return static_cast<std::string>(*this);
  }

  friend std::ostream& operator<<(std::ostream& os, const avx_bitset& bitset) {
    return os << static_cast<std::string>(bitset);
  }

  ~avx_bitset() { internal::Deallocate(data_); }

 private:
  avx_bitset(uint64_t* data, int size)
      : data_(data), n_(size), blocks_256_(size >> 8), blocks_64_(size >> 6) {}

  void MoveFrom(avx_bitset&& other) {
    data_ = std::exchange(other.data_, nullptr);
    n_ = std::exchange(other.n_, 0);
    blocks_256_ = std::exchange(other.blocks_256_, 0);
    blocks_64_ = std::exchange(other.blocks_64_, 0);
  }

  void CopyFrom(const avx_bitset& other) {
    data_ = internal::Allocate(other.n_);
    n_ = other.n_;
    blocks_256_ = other.blocks_256_;
    blocks_64_ = other.blocks_64_;
    for (int i = 0; i < internal::Ceil256(n_); ++i) {
      _mm256_stream_si256(
          reinterpret_cast<__m256i*>(data_ + (static_cast<ptrdiff_t>(i) * 4)),
          _mm256_stream_load_si256(reinterpret_cast<const __m256i*>(
              other.data_ + (static_cast<ptrdiff_t>(i) * 4))));
    }
  }

  void ShiftRightSmall(int shift) {
    if ((n_ & 63) != 0) {
      data_[blocks_64_] <<= shift;
      if (blocks_64_ == 0) {
        return;
      }
      data_[blocks_64_] |= data_[blocks_64_ - 1] >> (64 - shift);
    }
    for (int i = blocks_64_ - 1; i > 0; --i) {
      data_[i] <<= shift;
      data_[i] |= data_[i - 1] >> (64 - shift);
    }
    data_[0] <<= shift;
    RemoveTrash();
  }

  void ShiftLeftSmall(int shift) {
    for (int i = 0; i < blocks_64_ - 1; ++i) {
      data_[i] >>= shift;
      data_[i] |= data_[i + 1] << (64 - shift);
    }
    if (blocks_64_ != 0) {
      data_[blocks_64_ - 1] >>= shift;
      if ((n_ & 63) != 0) {
        data_[blocks_64_ - 1] |= data_[blocks_64_] << (64 - shift);
      }
    }
    if ((n_ & 63) != 0) {
      data_[blocks_64_] >>= shift;
    }
  }

  void RemoveTrash() {
    int id = internal::Blocks(n_) - 1;
    while ((id << 6) >= n_) {
      data_[id] = 0;
      id -= 1;
    }
    data_[id] &= internal::Mask(n_ & 63);
  }

  uint64_t* data_;

  int n_;
  int blocks_256_;
  int blocks_64_;
};

}  // namespace bitset
