#ifndef SEQUANT_EXTERNAL_INTERFACE_EXPORTSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_EXPORTSTEP_HPP

#include "execution_context.hpp"
#include "processing_step.hpp"

#include <SeQuant/core/space.hpp>

#include <nlohmann/json_fwd.hpp>

#include <filesystem>
#include <functional>
#include <map>
#include <string>
#include <string_view>
#include <variant>

namespace sequant::util::extint {

class ExportStep : public ProcessingStep {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

  std::size_t run(std::string_view step_id, ExecutionContext &ctx,
                  const std::vector<std::string_view> &inputs = {}) override;

 protected:
  struct ItfMeta {
    struct IndexSpaceMeta {
      std::string name;
      std::string tag;
    };

    std::map<IndexSpace, IndexSpaceMeta> index_spaces;
  };

  using meta_type = std::variant<std::monostate, ItfMeta>;

  std::string language_;
  bool optimize_ = true;
  std::filesystem::path filepath_;
  std::map<std::string, std::string, std::less<>> group_assoc_;
  meta_type meta_;
  std::vector<std::string> relative_order_;

  static meta_type parse_meta(std::string_view language,
                              const nlohmann::json &meta);

  friend class MetaAwareIftContext;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_EXPORTSTEP_HPP
