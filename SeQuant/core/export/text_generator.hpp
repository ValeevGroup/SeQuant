#ifndef SEQUANT_CORE_EXPORT_TEXTGENERATOR_HPP
#define SEQUANT_CORE_EXPORT_TEXTGENERATOR_HPP

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/generator.hpp>

namespace sequant {

struct TextGeneratorContext : ExportContext {};

class TextGenerator : public Generator<TextGeneratorContext> {
 public:
  TextGenerator() = default;
  ~TextGenerator();

  std::string get_format_name() const override;

  std::string represent(const Tensor &tensor,
                        const Context &ctx = {}) const override;
  std::string represent(const Variable &variable,
                        const Context &ctx = {}) const override;
  std::string represent(const Constant &constant,
                        const Context &ctx = {}) const override;

  void create(const Tensor &tensor, bool zero_init,
              const Context &ctx = {}) override;
  void load(const Tensor &tensor, bool set_to_zero,
            const Context &ctx = {}) override;
  void set_to_zero(const Tensor &tensor, const Context &ctx = {}) override;
  void unload(const Tensor &tensor, const Context &ctx = {}) override;
  void destroy(const Tensor &tensor, const Context &ctx = {}) override;
  void persist(const Tensor &tensor, const Context &ctx) override;

  void create(const Variable &variable, bool zero_init,
              const Context &ctx = {}) override;
  void load(const Variable &variable, bool set_to_zero,
            const Context &ctx = {}) override;
  void set_to_zero(const Variable &variable, const Context &ctx = {}) override;
  void unload(const Variable &variable, const Context &ctx = {}) override;
  void destroy(const Variable &variable, const Context &ctx = {}) override;
  void persist(const Variable &variable, const Context &ctx) override;

  void compute(const Expr &expression, const Variable &result,
               const Context &ctx = {}) override;
  void compute(const Expr &expression, const Tensor &result,
               const Context &ctx = {}) override;

  void declare(const Index &idx, const Context &ctx = {}) override;
  void declare(const Variable &variable, const Context &ctx) override;
  void declare(const Tensor &tensor, const Context &ctx = {}) override;

  void insert_comment(const std::string &comment,
                      const Context &ctx = {}) override;
  void insert_blank_lines(std::size_t count, const Context &ctx = {}) override;

  std::string get_generated_code() const override;

 private:
  std::string m_generated;

  std::string stringify(const Expr &expr) const;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_TEXTGENERATOR_HPP
