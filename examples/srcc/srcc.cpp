#include <SeQuant/core/op.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/formalism.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

#include <clocale>

using namespace sequant;
using namespace sequant::mbpt::sr;

namespace {

#define runtime_assert(tf)                                         \
  if (!(tf)) {                                                     \
    std::ostringstream oss;                                        \
    oss << "failed assert at line " << __LINE__ << " in function " \
        << __func__;                                               \
    throw std::runtime_error(oss.str().c_str());                   \
  }

TimerPool<32> tpool;

/// types of CC equations to solve
enum class EqnType { t, lambda };

/// maps equation type to string
inline const std::map<EqnType, std::wstring> type2str = {
    {EqnType::t, L"t"}, {EqnType::lambda, L"lambda"}};

/// maps equation type string to enum
inline const std::map<std::string, EqnType> str2type = {
    {"t", EqnType::t}, {"lambda", EqnType::lambda}};

/// maps unoccupied basis type string to enum
inline const std::map<std::string, mbpt::CSVFormalism> str2uocc = {
    {"std", mbpt::CSVFormalism::NonCSV}, {"csv", mbpt::CSVFormalism::CSV}};

// profiles evaluation of all CC equations for a given ex rank N with projection
// ex rank PMIN .. P
class compute_cceqvec {
  size_t P, PMIN, N;
  EqnType type;

 public:
  compute_cceqvec(size_t p, size_t pmin, size_t n, EqnType t = EqnType::t)
      : P(p), PMIN(pmin), N(n), type(t) {}

  void operator()(bool print, bool screen, bool use_topology,
                  bool use_connectivity, bool canonical_only) {
    tpool.start(N);
    std::vector<ExprPtr> eqvec;
    constexpr auto debug = true;
    switch (type) {
      case EqnType::t:
        eqvec = cceqs{N, P, PMIN}.t(screen, use_topology, use_connectivity,
                                    canonical_only);
        if constexpr (debug) {
          // troubleshoot against non-CSV CCSDTQ
          for (auto R = PMIN; R <= P; ++R) {
            std::wcout << "# of R" << R
                       << " terms in CSV-CC = " << eqvec.at(R)->size()
                       << std::endl;
          }
          if (N == 4) {
            auto r4_g_t3_t3_equiv_terms = std::make_shared<Sum>();
            ranges::for_each(eqvec.at(4)->expr(), [&](auto& summand) {
              if (summand->at(2).template as<Tensor>().rank() == 3 &&
                  summand->at(3).template as<Tensor>().rank() == 3 &&
                  summand.template as<Product>().scalar() == rational(1, 96)) {
                r4_g_t3_t3_equiv_terms->append(summand);
              }
            });
            std::wcout << "R4 CSV-CC terms that should have been reduced = "
                       << to_latex_align(r4_g_t3_t3_equiv_terms, 0, 1)
                       << std::endl;
            simplify(r4_g_t3_t3_equiv_terms);
            std::wcout << "simplify(r4_g_t3_t3_equiv_terms).size() = "
                       << r4_g_t3_t3_equiv_terms->size() << std::endl;

            // print automorphisms
            auto print_auts =
                [](const std::vector<std::vector<unsigned int>>& aut_generators,
                   auto&& stream, auto&& vlabels, bool use_labels) {
                  ranges::for_each(aut_generators, [&stream, &vlabels,
                                                    &use_labels](auto&& gen) {
                    // see bliss::print_permutation
                    auto print = [&stream, &vlabels, &use_labels](
                                     const std::vector<unsigned int>& perm) {
                      const unsigned int offset = 0;
                      const unsigned int N = perm.size();
                      for (unsigned int i = 0; i < N; i++) {
                        unsigned int j = perm[i];
                        if (j == i) continue;
                        bool is_first = true;
                        while (j != i) {
                          if (j < i) {
                            is_first = false;
                            break;
                          }
                          j = perm[j];
                        }
                        if (!is_first) continue;
                        stream << "("
                               << (use_labels ? vlabels.at(i)
                                              : std::to_wstring(i + offset))
                               << ",";
                        j = perm[i];
                        while (j != i) {
                          stream << (use_labels ? vlabels.at(j)
                                                : std::to_wstring(j + offset));
                          j = perm[j];
                          if (j != i) stream << ",";
                        }
                        stream << ")";
                      }
                    };

                    print(gen);
                    stream << std::endl;
                  });
                };

            for (auto& term : r4_g_t3_t3_equiv_terms->expr()) {
              std::wcout << "term = " << to_latex_align(term, 0, 1)
                         << std::endl;
              TensorNetwork tn(term->as<Product>().factors());
              auto [graph, vlabels, vcolors, vtypes] = tn.make_bliss_graph();

              std::basic_ostringstream<wchar_t> oss;
              graph->write_dot(oss, vlabels);
              std::wcout << "colored graph's dot = " << oss.str() << std::endl;

              bliss::Stats stats;
              graph->set_splitting_heuristic(bliss::Graph::shs_fsm);

              std::vector<std::vector<unsigned int>> aut_generators;
              auto save_aut = [&aut_generators](const unsigned int n,
                                                const unsigned int* aut) {
                aut_generators.emplace_back(aut, aut + n);
              };
              graph->find_automorphisms(
                  stats, &bliss::aut_hook<decltype(save_aut)>, &save_aut);
              std::basic_ostringstream<wchar_t> oss2;
              print_auts(aut_generators, oss2, vlabels, true);
              std::wcout << oss2.str() << std::endl;
            }
          }
          ranges::for_each(eqvec, [](ExprPtr& eq) {
            if (eq) {
              decsv(eq);
            }
          });
          auto resetter = set_scoped_default_formalism(
              mbpt::Formalism::make_default().set(mbpt::CSVFormalism::NonCSV));
          auto eqvec_ref = cceqs{N, P, PMIN}.t(
              screen, use_topology, use_connectivity, canonical_only);
          for (auto R = PMIN; R <= P; ++R) {
            std::wcout << "# of R" << R
                       << " terms in CC = " << eqvec_ref.at(R)->size()
                       << std::endl;
          }
          for (auto R = PMIN; R <= P; ++R) {
            auto diff = eqvec[R] - eqvec_ref[R];
            simplify(diff);
            std::wcout << "R" << R << "(expS" << N
                       << "): difference between standard and de-CSV'ed eqn "
                          "versions = "
                       << to_latex_align(diff, 0, 1) << std::endl;
            ranges::for_each(diff->expr(), [](const auto& summand) {
              std::wcout << "hash_value(" << to_latex(summand)
                         << ") = " << summand->hash_value() << std::endl;
              if (summand.template is<Product>()) {
                ranges::for_each(
                    summand.template as<Product>().factors(),
                    [](const auto& factor) {
                      std::wcout << "  hash_value(" << to_latex(factor)
                                 << ") = " << factor->hash_value() << std::endl;
                      if (factor.template is<Tensor>()) {
                        ranges::for_each(factor.template as<Tensor>().braket(),
                                         [](const auto& idx) {
                                           std::wcout
                                               << "    hash_value("
                                               << to_latex(idx)
                                               << ") = " << hash_value(idx)
                                               << std::endl;
                                         });
                      }
                    });
              }
            });
            // compare term by term
            std::wcout << "eqvec[R].size() = " << eqvec[R]->size() << std::endl;
            std::wcout << "eqvec_ref[R].size() = " << eqvec_ref[R]->size()
                       << std::endl;
            ranges::for_each(eqvec[R].as<Sum>().summands(), [&eqvec_ref,
                                                             R](auto&
                                                                    _summand) {
              auto summand = _summand->clone();
              // std::wcout << "summand = " << to_latex(summand) << std::endl;
              summand->rapid_canonicalize();
              canonicalize(summand);
              summand->rapid_canonicalize();
              // std::wcout << "canonicalize(summand) = " << to_latex(summand)
              // << std::endl;
              auto hash = summand->hash_value();
              auto it_by_hash = ranges::find_if(
                  eqvec_ref[R]->expr(),
                  [hash](const auto& s) { return s->hash_value() == hash; });
              auto it_by_hash_and_prefactor = ranges::find_if(
                  eqvec_ref[R]->expr(), [hash, &summand](const auto& s) {
                    return s->hash_value() == hash &&
                           s.template as<Product>().scalar() ==
                               summand.template as<Product>().scalar();
                  });
              if (it_by_hash == eqvec_ref[R]->expr().end()) {
                std::wcout << "summand " << to_latex(summand)
                           << " not found in reference eqn" << std::endl;
              } else if (it_by_hash_and_prefactor ==
                         eqvec_ref[R]->expr().end()) {
                std::wcout
                    << "summand " << to_latex(summand)
                    << " found in reference eqn with different prefactor, as "
                    << to_latex(*it_by_hash) << std::endl;
              } else {
                //                std::wcout << "summand " << to_latex(*summand)
                //                           << " found in reference eqn" <<
                //                           std::endl;
              }
            });
          }
        }
        break;
      case EqnType::lambda:
        eqvec = cceqs{N, P, PMIN}.lambda(screen, use_topology, use_connectivity,
                                         canonical_only);
        break;
    }
    tpool.stop(N);
    std::wcout << std::boolalpha << "CC equations [type=" << type2str.at(type)
               << ",rank=" << N << ",screen=" << screen
               << ",use_topology=" << use_topology
               << ",use_connectivity=" << use_connectivity
               << ",canonical_only=" << canonical_only << "] computed in "
               << tpool.read(N) << " seconds" << std::endl;
    for (size_t R = PMIN; R <= P; ++R) {
      std::wcout << "R" << R << "(expS" << N << ") has " << eqvec[R]->size()
                 << " terms:" << std::endl;
      if (print) std::wcout << to_latex_align(eqvec[R], 20, 3) << std::endl;

      // validate known sizes of some CC residuals
      // N.B. # of equations depends on whether we use symmetric or
      // antisymmetric amplitudes
      if (mbpt::get_default_formalism().two_body_interaction() ==
          mbpt::TwoBodyInteraction::Antisymm) {
        if (type == EqnType::t) {
          if (R == 1 && N == 1) runtime_assert(eqvec[R]->size() == 8);
          if (R == 1 && N == 2) runtime_assert(eqvec[R]->size() == 14);
          if (R == 2 && N == 2) runtime_assert(eqvec[R]->size() == 31);
          if (R == 1 && N == 3) runtime_assert(eqvec[R]->size() == 15);
          if (R == 2 && N == 3) runtime_assert(eqvec[R]->size() == 37);
          if (R == 3 && N == 3) runtime_assert(eqvec[R]->size() == 47);
          if (R == 4 && N == 4) runtime_assert(eqvec[R]->size() == 74);
          if (R == 5 && N == 5) runtime_assert(eqvec[R]->size() == 99);
        }
      } else {
        if (type == EqnType::t) {
          if (R == 1 && N == 2) runtime_assert(eqvec[R]->size() == 26);
          if (R == 2 && N == 2) runtime_assert(eqvec[R]->size() == 55);
          if (R == 1 && N == 3) runtime_assert(eqvec[R]->size() == 30);
          if (R == 2 && N == 3) runtime_assert(eqvec[R]->size() == 73);
          if (R == 3 && N == 3) runtime_assert(eqvec[R]->size() == 93);
          if (R == 4 && N == 4) runtime_assert(eqvec[R]->size() == 149);
        }
      }
    }
  }

 private:
  static void decsv(ExprPtr& eq) {
    auto decsv_visitor = [](ExprPtr& expr) -> void {
      if (expr.is<Tensor>()) {  // remove protoindices
        const auto& tensor = expr.as<Tensor>();
        // visitation is depth-first, must skip overlaps
        if (ranges::any_of(tensor.braket(),
                           [](auto& idx) { return idx.has_proto_indices(); })) {
          auto protoindexfree_bra = tensor.bra() |
                                    ranges::views::transform([](auto& idx) {
                                      return idx.drop_proto_indices();
                                    }) |
                                    ranges::to<container::svector<Index>>;
          auto protoindexfree_ket = tensor.ket() |
                                    ranges::views::transform([](auto& idx) {
                                      return idx.drop_proto_indices();
                                    }) |
                                    ranges::to<container::svector<Index>>;
          auto protoindexfree_tensor = std::make_shared<Tensor>(
              tensor.label(), protoindexfree_bra, protoindexfree_ket,
              tensor.symmetry(), tensor.braket_symmetry(),
              tensor.particle_symmetry());
          expr = protoindexfree_tensor;
        }
      }
    };
    eq->visit(decsv_visitor, /* atoms_only = */ false);
    // std::wcout << "decsv(eq) = " << to_latex(eq) << std::endl;
    FWickTheorem wk(eq);
    wk.reduce(eq);  // this reduces the overlaps
    std::wcout << "decsv+reduce(eq).size() = " << eq->size() << std::endl;
    //    simplify(eq);
    //    std::wcout << "decsv+reduce+simplify(eq).size() = " << eq->size() <<
    //    std::endl;
  };

};  // class compute_cceqvec

// profiles evaluation of all CC equations with ex rank 2 .. N
class compute_all {
  size_t NMAX;
  EqnType type;

 public:
  compute_all(size_t nmax, EqnType t = EqnType::t) : NMAX(nmax), type(t) {}

  void operator()(bool print = true, bool screen = true,
                  bool use_topology = true, bool use_connectivity = true,
                  bool canonical_only = true) {
    for (size_t N = 1; N <= NMAX; ++N)
      compute_cceqvec{N, 1, N, type}(print, screen, use_topology,
                                     use_connectivity, canonical_only);
  }
};  // class compute_all

}  // namespace

int main(int argc, char* argv[]) {
  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  std::wcout.sync_with_stdio(false);
  std::wcerr.sync_with_stdio(false);
  std::wcout.imbue(std::locale("en_US.UTF-8"));
  std::wcerr.imbue(std::locale("en_US.UTF-8"));
  std::wcout.sync_with_stdio(true);
  std::wcerr.sync_with_stdio(true);
  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::set_default_context(
      SeQuant(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
              BraKetSymmetry::conjugate, SPBasis::spinorbital));
  mbpt::set_default_convention();
  mbpt::set_default_formalism();

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  // set_num_threads(1);

#ifndef NDEBUG
  const size_t DEFAULT_NMAX = 3;
#else
  const size_t DEFAULT_NMAX = 4;
#endif
  const size_t NMAX = argc > 1 ? std::atoi(argv[1]) : DEFAULT_NMAX;
  const std::string eqn_type_str = argc > 2 ? argv[2] : "t";
  const EqnType eqn_type = str2type.at(eqn_type_str);
  const std::string uocc_type_str = argc > 3 ? argv[3] : "std";
  const mbpt::CSVFormalism uocc_type = str2uocc.at(uocc_type_str);
  auto resetter = set_scoped_default_formalism(
      mbpt::Formalism::make_default().set(uocc_type));

  // change to true to print out the resulting equations
  constexpr bool print = false;
  // change to true to print stats
  Logger::get_instance().wick_stats = false;

  ranges::for_each(std::array<bool, 2>{false, true}, [=](const bool screen) {
    ranges::for_each(
        std::array<bool, 2>{false, true}, [=](const bool use_topology) {
          ranges::for_each(
              std::array<bool, 2>{false, true}, [=](const bool canonical_only) {
                tpool.clear();
                // comment out to run all possible combinations
                if (screen && use_topology && canonical_only)
                  compute_all{NMAX, eqn_type}(print, screen, use_topology, true,
                                              canonical_only);
              });
        });
  });
}
