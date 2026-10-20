//
// Plain data types of the JSON "schedule document": the reproducible record of
// a batched-evaluation schedule (fingerprint, per-term tree shapes and node
// annotations, the value/block DAG and an optional run section), plus its
// reader, writer, shape validation and a structural comparison.
//

#ifndef SEQUANT_EVAL_SCHEDULE_DOCUMENT_HPP
#define SEQUANT_EVAL_SCHEDULE_DOCUMENT_HPP

#include <cstddef>
#include <iosfwd>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace sequant {
class Index;
}  // namespace sequant

namespace sequant::eval::schedule {

/// the spelling of \p index in every axis label and mode list of a schedule
/// document: the V1 serializer's (proto indices joined by ","), as it appears
/// inside the term's factor strings
std::string index_label(Index const& index);

/// value of the `format` key of every schedule document
inline constexpr char const* kFormat = "sequant.schedule/1";

/// mirrors \c BatchModeType
enum class ModeKind { Contracted, External };

/// a batched axis; label = index_label(axis index)
struct AxisRecord {
  std::string label;
  ModeKind kind;
};

/// per-contraction-node batch annotation
struct NodeAnnotationRecord {
  std::vector<AxisRecord> axes, opened_here;
  bool order_aware = false;
};

/// binary tree shape: a leaf carries an ordinal (into the term's factors), an
/// inner node has two children. Deep-copyable.
struct TreeShape {
  std::optional<std::size_t> leaf;
  std::unique_ptr<TreeShape> left, right;

  TreeShape() = default;
  TreeShape(std::size_t ordinal) : leaf(ordinal) {}
  TreeShape(TreeShape const& other)
      : leaf(other.leaf),
        left(other.left ? std::make_unique<TreeShape>(*other.left) : nullptr),
        right(other.right ? std::make_unique<TreeShape>(*other.right)
                          : nullptr) {}
  TreeShape(TreeShape&&) noexcept = default;
  TreeShape& operator=(TreeShape const& other) {
    if (this != &other) {
      TreeShape tmp(other);
      *this = std::move(tmp);
    }
    return *this;
  }
  TreeShape& operator=(TreeShape&&) noexcept = default;
};

struct TermRecord {
  std::string key;
  std::vector<std::string> factors;
  std::string scalar;
  std::optional<TreeShape> tree;  // nullopt = pure scalar term
  std::vector<NodeAnnotationRecord> nodes;
};

struct RootRecord {
  std::string label;
  std::vector<TermRecord> terms;
};

struct ContextRecord {
  std::string vacuum, spbasis;
  bool strict_braket_symmetry = false;
  std::vector<std::string> cardinal_tensor_labels, spaces;
};

struct Fingerprint {
  std::vector<std::string> roots;
  std::map<std::string, std::pair<bool, bool>>
      batchable;  // space -> (contracted, external)
  std::vector<std::string> volatile_leaves;
  bool persistent_only = false, batch_spectator_indices = false;
  double accumulation_factor = 1;
  std::string scheduler, objective_function, reorder;
  std::map<std::string, std::string> cse;  // CSEOptions fields as strings
  std::vector<std::string> mode_order;
  ContextRecord context;
};

struct OccurrenceRecord {
  std::string term;
  std::size_t node = 0;
  std::optional<std::size_t> consumer;
  std::vector<std::string> home, enclosing;
  std::vector<AxisRecord> opens;
};

/// role: LoopLocal|Reduction|LoopCarried|LoopInvariant
struct PerAxisRecord {
  std::string mode;
  std::string role;
};

struct ValueRecord {
  std::string id, expr;
  std::string hash, key;
  bool leaf = false;
  std::vector<std::string> result_modes, home_modes, enclosing_modes,
      divergent_modes;
  std::size_t home_block = 0;
  int home_depth = 0;
  bool persistent = false;
  std::size_t first_use = 0, last_use = 0;
  std::vector<std::string> build_site, home_floor, forced_split_axes;
  std::vector<PerAxisRecord> per_axis;
  std::vector<OccurrenceRecord> occurrences;
  std::vector<std::string> siblings;
};

struct StepRecord {
  std::optional<std::string> build;
  std::optional<std::size_t> block;
};

struct BlockRecord {
  std::size_t id = 0;
  std::optional<std::size_t> parent;
  std::optional<std::string> axis;
  std::string kind;  // "root" | "contracted" | "external"
  std::size_t depth = 0;
  int slot = 0, altitude = 0, pass = 0;
  std::vector<StepRecord> steps;
  /// (value id, "AccumulateSum"|"AccumulateScatter")
  std::vector<std::pair<std::string, std::string>> outputs;
};

struct LoopRecord {
  std::string space;
  std::size_t depth = 0;
  std::string kind;
};

struct EdgeRecord {
  std::string from, to;
  bool sliced = false;
};

struct DagRecord {
  std::vector<ValueRecord> values;
  std::vector<BlockRecord> blocks;
  std::vector<LoopRecord> loops;
  std::vector<std::tuple<LoopRecord, LoopRecord, std::size_t>> loop_order;
  std::vector<EdgeRecord> edges;
};

struct RunRecord {
  std::map<std::string, std::size_t> extents, batch_targets;  // by space
  std::map<std::size_t, std::size_t> n_batches;  // block id -> count
  /// value id -> bytes (largest result a production left in its cell)
  std::map<std::string, std::size_t> footprint;
  std::map<std::string, std::size_t> builds;  // value id -> build count
  /// optimizer settings excluded from the fingerprint
  std::map<std::string, std::string> planner;
};

struct Document {
  std::string producer_version, producer_date;
  Fingerprint fingerprint;
  std::vector<RootRecord> roots;
  std::optional<DagRecord> dag;
  std::optional<RunRecord> run;
};

/// pretty (2-space indent), keys in struct order
std::string to_json(Document const&);
void write(Document const&, std::ostream&);

/// \throws sequant::Exception on malformed or invalid documents
Document from_json(std::string_view);
Document read(std::istream&);

struct Difference {
  std::string path, left, right;
};

/// structural comparison; ignores `run`, `producer_*` and, on DAG values,
/// `hash`, `key`, `first_use`, `last_use`
std::vector<Difference> compare(Document const& a, Document const& b);

/// "roots=<n> terms=<n> values=<n> blocks=<n> depth=<n>"; the last three are
/// "-" when the document has no `dag` section. depth = the deepest nesting of
/// the block tree (number of ancestors of a block), not BlockRecord::depth
std::string summary_line(Document const&);

}  // namespace sequant::eval::schedule

#endif  // SEQUANT_EVAL_SCHEDULE_DOCUMENT_HPP
