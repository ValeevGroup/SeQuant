#ifndef SEQUANT_EXTERNAL_INTERFACE_PROCESSINGSTEP_HPP
#define SEQUANT_EXTERNAL_INTERFACE_PROCESSINGSTEP_HPP

#include <SeQuant/core/expr.hpp>

#include <nlohmann/json_fwd.hpp>

#include <memory>
#include <string>
#include <vector>
#include <filesystem>


namespace sequant::util::extint {

struct ProcessingInput {
	virtual ~ProcessingInput() = default;
};

struct VoidInput : public ProcessingInput {};

struct ExpressionInput : public ProcessingInput {
	std::vector<ResultExpr> expressions;
};

struct ExportTreeInput : public ProcessingInput {
};

struct ProcessingOutput {
	virtual ~ProcessingOutput() = default;
};

struct VoidOutput : public ProcessingOutput {};

struct ExpressionOutput : public ProcessingOutput {
	std::vector<ResultExpr> expressions;
};

struct ExportTreeOutput : public ProcessingOutput {
};


class ProcessingStep {
public:
	ProcessingStep() = default;
	virtual ~ProcessingStep() = default;

	virtual std::string kind() const = 0;

	virtual bool accepts_options() const = 0;
	virtual bool requires_options() const = 0;
	virtual void set_options(const nlohmann::json &options) = 0;

	virtual std::unique_ptr<ProcessingOutput> run(const ProcessingInput &input); 
};


class ReadInputStep : public ProcessingStep {
public:
	std::string kind() const override;

	bool accepts_options() const override;
	bool requires_options() const override;
	void set_options(const nlohmann::json &options) override;

	std::unique_ptr<ProcessingOutput> run(const ProcessingInput &) override;

private:
	std::vector<std::filesystem::path> input_paths_;
};


}

#endif // SEQUANT_EXTERNAL_INTERFACE_PROCESSINGSTEP_HPP
