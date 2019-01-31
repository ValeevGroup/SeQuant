//
// Created by Eduard Valeyev on 2019-01-30.
//

#ifndef SEQUANT2_TAG_HPP
#define SEQUANT2_TAG_HPP

#include <boost/any.hpp>

namespace sequant2 {

class Taggable {
 public:
  using any = boost::any;

  Taggable() noexcept : tag_{} {
    assert(!has_tag());
  }

  /// tags this object with tag @c t
  template<typename T>
  void tag(const T &t) {
    assert(tag_.empty());
    tag_ = t;
    assert(!tag_.empty());
  }

  /// returns this object's tag
  template<typename T>
  const T &tag() const {
    assert(!tag_.empty());
    using boost::any_cast;
    return any_cast<T>(tag_);
  }

  /// @return true if tag has been assigned
  bool has_tag() const { return !tag_.empty(); }

  /// resets this object's tag
  void reset_tag() { tag_ = any{}; }

 private:
  any tag_;
};

}  // namespace sequant2

#endif //SEQUANT2_TAG_HPP
