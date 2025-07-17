#ifndef SEQUANT_EXPRESSIONS_LABELED_HPP
#define SEQUANT_EXPRESSIONS_LABELED_HPP

#include <string_view>

namespace sequant {

/// An object with a string label that be used for defining a canonical order of
/// expressions (defined at runtime). These labels are not in general
/// object-unique, i.e. 2 objects that are not equal can return same label.
class Labeled {
 public:
  Labeled() = default;
  virtual ~Labeled() = default;

  /// @return object's label. 2 objects that are not equal can return same
  /// label, i.e. the label does not identify the object uniquely. Thus labels
  /// can be viewed as immutable string tags that for some object types can be
  /// used for sorting and other purposes.
  /// @sa to_latex() for producing unique string representation
  virtual std::wstring_view label() const = 0;
};

/// A Labeled object that also allows changing the label
class MutatableLabeled : public Labeled {
 public:
  using Labeled::Labeled;

  /// Sets the label of this object
  /// @param label The new label of this object
  virtual void set_label(std::wstring label) = 0;
};

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_LABELED_HPP
