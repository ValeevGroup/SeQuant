#include <SeQuant/core/bliss.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/tensor_network_v2.hpp>
#include <SeQuant/core/tensor_network_v3.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <cassert>
#include <iostream>
#include <optional>
#include <string>

using namespace sequant;

template <typename TN>
std::optional<TN> make_tn(const ExprPtr &expr) {
  if (expr.is<Tensor>()) {
    return TN({expr});
  } else if (expr.is<Product>()) {
    for (const ExprPtr &factor : expr.as<Product>().factors()) {
      if (!factor.is<Tensor>()) {
        return {};
      }
    }

    return TN(expr.as<Product>().factors());
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
  std::wcout
      << "  --v1       Use original TensorNetwork (aka TensorNetworkV1)\n";
  std::wcout << "  --v2       Use TensorNetworkV2\n";
  std::wcout << "  --v3       Use TensorNetworkV3 [default]\n";
  std::wcout << "  --no-named Treat all indices as unnamed (even if they are "
                "external)\n";
}

int main(int argc, char **argv) {
  set_locale();
  sequant::set_default_context(
      {.index_space_registry_shared_ptr = mbpt::make_sr_spaces(),
       .vacuum = Vacuum::SingleProduct});

  bool use_named_indices = true;
  int version = 3;
  const TensorNetwork::named_indices_t empty_named_indices;

  if (argc <= 1) {
    print_help();
    return 0;
  }

  for (std::size_t i = 1; i < static_cast<std::size_t>(argc); ++i) {
    std::wstring current = toUtf16(argv[i]);
    if (current == L"--help") {
      print_help();
      return 0;
    } else if (current == L"--no-named") {
      use_named_indices = false;
      continue;
    } else if (current == L"--v1") {
      version = 1;
      continue;
    } else if (current == L"--v2") {
      version = 2;
      continue;
    } else if (current == L"--v3") {
      version = 3;
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

    if (version == 1) {
      std::optional<TensorNetwork> network = make_tn<TensorNetwork>(expr);
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
      graph->write_dot(std::wcout, {.labels = vlabels});
    } else {
      auto make_graph = [&](auto *tn_ptr) -> int {
        using TN = std::decay_t<decltype(*tn_ptr)>;
        std::optional<TN> network = make_tn<TN>(expr);
        if (!network.has_value()) {
          std::wcout << "Failed to construct tensor network for input '"
                     << to_latex(expr) << "'" << std::endl;
          return 3;
        }

        auto graph = network->create_graph(
            {.named_indices =
                 use_named_indices ? nullptr : &empty_named_indices});
        std::wcout << "Graph for '" << to_latex(expr) << "'\n";
        graph.bliss_graph->write_dot(std::wcout,
                                     {.labels = graph.vertex_labels});
        return 0;
      };
      switch (version) {
        case 2:
          return make_graph(static_cast<TensorNetworkV2 *>(nullptr));
        case 3:
          return make_graph(static_cast<TensorNetworkV3 *>(nullptr));
        default:
          abort();  // unreachable
      }
    }
  }
}
