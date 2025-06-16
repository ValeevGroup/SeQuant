#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/tensor_network_v2.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <cassert>
#include <codecvt>
#include <iostream>
#include <locale>
#include <optional>
#include <string>

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

std::optional<TensorNetworkV2> to_network_v2(const ExprPtr &expr) {
  if (expr.is<Tensor>()) {
    return TensorNetworkV2({expr});
  } else if (expr.is<Product>()) {
    for (const ExprPtr &factor : expr.as<Product>().factors()) {
      if (!factor.is<Tensor>()) {
        return {};
      }
    }

    return TensorNetworkV2(expr.as<Product>().factors());
  } else {
    return {};
  }
}

void print_help() {
  std::wcout << "Helper to generate dot (GraphViz) representations of tensor "
                "network graphs.\n";
  std::wcout << "Usage:\n";
  std::wcout
      << "  <exe> [options] <network 1> [<network 2> [... [<network N>] ] ]\n";
  std::wcout << "Options:\n";
  std::wcout << "  --help     Shows this help message\n";
  std::wcout << "  --v2       Use TensorNetworkV2\n";
  std::wcout << "  --no-named Treat all indices as unnamed (even if they are "
                "external)\n";
}

int main(int argc, char **argv) {
  set_locale();
  sequant::set_default_context(Context(
      mbpt::make_sr_spaces(), Vacuum::SingleProduct, IndexSpaceMetric::Unit,
      BraKetSymmetry::conjugate, SPBasis::spinorbital));

  bool use_named_indices = true;
  bool use_tnv2 = false;
  const TensorNetwork::named_indices_t empty_named_indices;

  if (argc <= 1) {
    print_help();
    return 0;
  }

  for (std::size_t i = 1; i < static_cast<std::size_t>(argc); ++i) {
    std::wstring current = from_utf8(argv[i]);
    if (current == L"--help") {
      print_help();
      return 0;
    } else if (current == L"--no-named") {
      use_named_indices = false;
      continue;
    } else if (current == L"--v2") {
      use_tnv2 = true;
      continue;
    }

    ExprPtr expr;
    try {
      expr = parse_expr(current);
    } catch (const ParseError &e) {
      std::wcout << "Failed to parse expression '" << current
                 << "': " << e.what() << std::endl;
      return 1;
    }
    assert(expr);

    if (!use_tnv2) {
      std::optional<TensorNetwork> network = to_network(expr);
      if (!network.has_value()) {
        std::wcout << "Failed to construct tensor network for input '"
                   << to_latex(expr) << "'" << std::endl;
        return 2;
      }

      auto [graph, vlabels, vtexlabels, vcolors, vtypes] =
          network->make_bliss_graph(
              {.named_indices =
                   (use_named_indices ? nullptr : &empty_named_indices)});
      std::wcout << "Graph for '" << to_latex(expr) << "'\n";
      graph->write_dot(std::wcout, vlabels);
    } else {
      std::optional<TensorNetworkV2> network = to_network_v2(expr);
      if (!network.has_value()) {
        std::wcout << "Failed to construct tensor network for input '"
                   << to_latex(expr) << "'" << std::endl;
        return 2;
      }

      TensorNetworkV2::Graph graph = network->create_graph(
          {.named_indices =
               use_named_indices ? nullptr : &empty_named_indices});
      std::wcout << "Graph for '" << to_latex(expr) << "'\n";
      graph.bliss_graph->write_dot(std::wcout, graph.vertex_labels);
    }
  }
}
