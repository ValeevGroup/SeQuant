#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <string>
#include <cassert>
#include <codecvt>
#include <locale>
#include <iostream>
#include <optional>

	using namespace sequant;

std::wstring from_utf8(std::string_view str) {
	std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
	return converter.from_bytes(std::string(str));
}

std::optional<TensorNetwork> to_network(const ExprPtr &expr) {
	if (expr.is<Tensor>()) {
		return TensorNetwork({expr});
	} else if (expr.is<Product>()) {
		for (const ExprPtr &factor : expr.as<Product>().factors()) {
			if (!factor.is<Tensor>()) {
				return {};
			}
		}

		return TensorNetwork(expr.as<Product>().factors());
	} else {
		return {};
	}
}

int main(int argc, char **argv) {

	set_locale();
	mbpt::set_default_convention();

	for (std::size_t i = 1; i < static_cast<std::size_t>(argc); ++i) {
		std::wstring current = from_utf8(argv[i]);
		ExprPtr expr;
		try {
			expr = parse_expr(current);
		} catch (const ParseError &e) {
			std::wcout << "Failed to parse expression '" << current << "': " << e.what() << std::endl;
			return 1;
		}
		assert(expr);

		std::optional<TensorNetwork> network = to_network(expr);
		if (!network.has_value()) {
			std::wcout << "Failed to construct tensor network for input '" << current << "'" << std::endl;
			return 2;
		}

		TensorNetwork::Graph graph = network->create_graph();
		std::wcout << "Graph for '" << current << "'\n";
		graph.bliss_graph->write_dot(std::wcout, graph.vertex_labels);
	}
}
