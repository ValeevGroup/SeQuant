#ifndef SEQUANT_EXTERNAL_INTERFACE_FILTERSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_FILTERSTEP_HPP

#include "execution_context.hpp"
#include "processing_data.hpp"
#include "processing_step.hpp"

#include <SeQuant/core/expr_fwd.hpp>

#include <nlohmann/json_fwd.hpp>

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

class ExpressionFilter {
 public:
  struct Rule {
    bool negate = false;

    virtual bool matches(const Expr &expr) const = 0;
  };

  ExpressionFilter() = default;

  bool matches(const Expr &) const;

  void set_require_all(bool require);
  void add_rule(std::unique_ptr<Rule> rule);

 private:
  bool require_all_ = true;
  std::vector<std::unique_ptr<Rule>> rules_;
};

class FilterStep : public OneByOneProcessingStep<ExpressionData, false> {
 public:
  std::string kind() const override;

  bool accepts_options() const override;
  bool requires_options() const override;
  void set_options(const nlohmann::json &options) override;

 protected:
  std::size_t process(std::string_view id_prefix, std::size_t id_start,
                      ExecutionContext &ctx,
                      const ExpressionData &data) override;

 private:
  std::map<std::string, ExpressionFilter, std::less<>> groups_;
  bool keep_empty_ = true;
};

}  // namespace sequant::util::extint

#endif  // SEQUANT_EXTERNAL_INTERFACE_FILTERSTEP_HPP
