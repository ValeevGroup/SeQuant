#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/eqs/cceqs.hpp>
#include <SeQuant/domain/eval/eval.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/optimize/optimize.hpp>
#include <SeQuant/domain/utils/binary_node.hpp>
#include <SeQuant/domain/utils/eval_expr.hpp>
#include <SeQuant/domain/utils/parse_expr.hpp>

#include <btas/btas.h>
#include <btas/tensor_func.h>
#include <fstream>
#include <range/v3/all.hpp>

#include <chrono>
#include <iostream>

auto const assert_imaginary_zero(sequant::Constant const& c) {
  assert(c.value().imag() == 0 && "complex scalar unsupported for real tensor");
}

std::array<size_t, 3> read_header(std::ifstream& fs) {
  size_t rank = 0, nocc = 0, nvirt = 0;
  std::string header{};
  std::getline(fs, header);
  auto hs = std::istringstream{header};
  hs >> rank;
  hs >> nocc;
  hs >> nvirt;
  return {rank, nocc, nvirt};
}  // read_header

auto read_header(std::string_view fname) {
  auto ifs = std::ifstream{fname};
  return read_header(ifs);
}  // read_header

auto read_tensor(std::string_view fname) {
  auto ifs = std::ifstream{fname};
  auto const [rank, nocc, nvirt] = read_header(ifs);
  auto range = btas::Range{ranges::views::repeat(nocc + nvirt) |
                           ranges::views::take(rank) | ranges::to_vector};
  auto tensor = btas::Tensor<double>{range};
  tensor.generate([&ifs]() {
    double x;
    ifs >> x;
    return x;
  });
  return tensor;
}  // read_tensor

auto const norm = [](btas::Tensor<double> const& tensor) {
  return std::sqrt(btas::dot(tensor, tensor));
};

template <typename Tensor_t>
class yield_leaf {
 private:
  Tensor_t const &fock, &eri, &t_vo, &t_vvoo;

  size_t const nocc, nvirt;

  auto range1_limits(sequant::Tensor const& tensor) {
    return tensor.const_braket() |
           ranges::views::transform([this](auto const& idx) {
             auto ao = sequant::IndexSpace::active_occupied;
             auto au = sequant::IndexSpace::active_unoccupied;
             auto sp = idx.space();
             assert(sp == ao || sp == au);

             return sp == ao ? nocc : nvirt;
           });
  }

 public:
  yield_leaf(size_t no, size_t nv, Tensor_t const& F, Tensor_t const& G,
             Tensor_t const& ampl_vo, Tensor_t const& ampl_vvoo)
      : nocc{no},
        nvirt{nv},
        fock{F},
        eri{G},
        t_vo{ampl_vo},
        t_vvoo{ampl_vvoo} {}

  Tensor_t operator()(sequant::Tensor const& texpr) {
    auto rank = texpr.bra_rank() + texpr.ket_rank();

    if (texpr.label() == L"t") {
      assert(rank == 2 || rank == 4 && "only t_vo and t_vvoo supported");
      return rank == 2 ? t_vo : t_vvoo;
    }

    assert((texpr.label() == L"g" || texpr.label() == L"f") &&
           "unsupported tensor label encountered");

    auto&& big_tensor = texpr.label() == L"g" ? eri : fock;

    auto r1_limits = range1_limits(texpr);
    auto iter_limits = r1_limits | ranges::views::transform([this](auto x) {
                         return x == nocc ? std::pair{size_t{0}, nocc}
                                          : std::pair{nocc, nocc + nvirt};
                       });

    auto slice = Tensor_t{btas::Range{r1_limits | ranges::to_vector}};

    if (iter_limits.size() == 2) {
      auto loop1 = iter_limits[0];
      auto loop2 = iter_limits[1];
      for (auto i = loop1.first; i < loop1.second; ++i)
        for (auto j = loop2.first; j < loop2.second; ++j)
          slice(i - loop1.first, j - loop2.first) = big_tensor(i, j);

    } else {  // iter_limits.size() == 4 true
      auto loop1 = iter_limits[0];
      auto loop2 = iter_limits[1];
      auto loop3 = iter_limits[2];
      auto loop4 = iter_limits[3];
      for (auto i = loop1.first; i < loop1.second; ++i)
        for (auto j = loop2.first; j < loop2.second; ++j)
          for (auto k = loop3.first; k < loop3.second; ++k)
            for (auto l = loop4.first; l < loop4.second; ++l)
              slice(i - loop1.first,    //
                    j - loop2.first,    //
                    k - loop3.first,    //
                    l - loop4.first) =  //
                  big_tensor(i, j, k, l);
    }

    return slice;
  }
};  // tensor yield

template <typename Tensor_t>
Tensor_t inode_evaluate_btas(
    sequant::utils::binary_node<sequant::utils::eval_expr> const& node,
    Tensor_t const& leval, Tensor_t const& reval) {
  assert((node->op() == sequant::utils::eval_expr::eval_op::Sum ||
          node->op() == sequant::utils::eval_expr::eval_op::Prod) &&
         "unsupported intermediate operation");

  assert_imaginary_zero(node.left()->scalar());
  assert_imaginary_zero(node.right()->scalar());

  auto index_hash = [](auto const& bk) {
    return ranges::views::transform(bk,
                                    [](auto const& idx) {
                                      size_t seed = 0;
                                      sequant::hash::combine(seed, idx.label());
                                      return static_cast<long>(seed);  //
                                    }) |
           ranges::to_vector;
  };  // index_hash

  if (node->op() == sequant::utils::eval_expr::eval_op::Prod) {
    auto scalar = node.left()->scalar().value().real() *
                  node.right()->scalar().value().real();

    auto lannot = index_hash(node.left()->tensor().const_braket());
    auto rannot = index_hash(node.right()->tensor().const_braket());
    auto this_annot = index_hash(node->tensor().const_braket());

    Tensor_t prod;

    btas::contract(     //
        scalar,         //
        leval, lannot,  //
        reval, rannot,  //
        0.0,            //
        prod, this_annot);
    return prod;
  }

  else {  // sum

    auto const post_annot = index_hash(node->tensor().const_braket());
    auto permute_and_scale = [&post_annot, &index_hash](auto const& btensor,
                                                        auto const& child_seqt,
                                                        auto scal) {
      auto pre_annot = index_hash(child_seqt.const_braket());
      Tensor_t result;
      btas::permute(btensor, pre_annot, result, post_annot);
      btas::scal(scal, result);
      return result;
    };

    auto lscal = node.left()->scalar().value().real();
    auto rscal = node.right()->scalar().value().real();

    auto sum = permute_and_scale(leval, node.left()->tensor(), lscal);
    sum += permute_and_scale(reval, node.right()->tensor(), rscal);

    return sum;
  }
}

template <typename Tensor_t>
Tensor_t evaluate_btas(
    sequant::utils::binary_node<sequant::utils::eval_expr> const& node,
    yield_leaf<Tensor_t>& yielder,
    sequant::eval::cache_manager<Tensor_t>& cman) {
  if (auto const& exists = cman.access(node->hash()); exists)
    return exists.value();
  if (node.leaf() && (node->tensor().label() != L"I"))
    return cman.store(node->hash(), yielder(node->tensor()));
  assert((!node.leaf()) && "this shouldn't happen");
  auto eval_result =
      inode_evaluate_btas(node, evaluate_btas(node.left(), yielder, cman),
                          evaluate_btas(node.right(), yielder, cman));
  return cman.store(node->hash(), std::move(eval_result));
}

struct eval_instance {
  sequant::utils::binary_node<sequant::utils::eval_expr> const& node;

  template <typename Tensor_t>
  auto evaluate(yield_leaf<Tensor_t>& yielder,
                sequant::eval::cache_manager<Tensor_t>& cman) {
    auto result = evaluate_btas(node, yielder, cman);

    btas::scal(node->scalar().value().real(), result);

    return result;
  }

  template <typename Tensor_t>
  auto evaluate_symm(yield_leaf<Tensor_t>& yielder,
                     sequant::eval::cache_manager<Tensor_t>& cman) {
    auto pre_symm = evaluate(yielder, cman);
    auto result = decltype(pre_symm){pre_symm.range()};
    result.fill(0.);

    auto const lannot = ranges::views::iota(size_t{0}, pre_symm.rank()) |
                        ranges::to<sequant::eval::perm_type>;

    auto symm_impl = [&result, &pre_symm, &lannot](auto const& annot) {
      decltype(result) temp;
      btas::permute(pre_symm, lannot, temp, annot);
      result += temp;
    };

    sequant::eval::symmetrize_tensor(pre_symm.rank(), symm_impl);
    return result;
  }

  template <typename Tensor_t>
  auto evaluate_asymm(yield_leaf<Tensor_t>& yielder,
                      sequant::eval::cache_manager<Tensor_t>& cman) {
    auto pre_asymm = evaluate(yielder, cman);
    auto result = decltype(pre_asymm){pre_asymm.range()};
    result.fill(0.);

    auto const lannot = ranges::views::iota(size_t{0}, pre_asymm.rank()) |
                        ranges::to<sequant::eval::perm_type>;

    auto asymm_impl = [&result, &pre_asymm,
                       &lannot](auto const& pwp) {  // pwp = phase with perm
      decltype(result) temp;
      btas::permute(pre_asymm, lannot, temp, pwp.perm);
      btas::scal(pwp.phase, temp);
      result += temp;
    };

    sequant::eval::antisymmetrize_tensor(pre_asymm.rank(), asymm_impl);
    return result;
  }

};  // evaluate_btas

int main(int argc, char** argv) {
  using std::cout;
  using std::endl;

  auto fock_input = "scratch/fock_h2o.dat";
  auto eri_input = "scratch/eri_h2o.dat";

  auto fock_head = read_header(fock_input);
  auto eri_head = read_header(eri_input);

  assert(fock_head[0] == 2 && "Fock tensor rank != 2");
  assert(eri_head[0] == 4 && "ERI tensor rank != 4");

  for (auto&& pair :
       ranges::views::zip(fock_head, eri_head) | ranges::views::tail)
    assert(pair.first == pair.second &&
           "incompatible Range1 objects in read tensors");

  auto const fock = read_tensor(fock_input);
  auto const eri = read_tensor(eri_input);

  size_t const nocc = fock_head[1];
  size_t const nvirt = fock_head[2];

  auto t_vo = btas::Tensor<double>{btas::Range{nvirt, nocc}};
  auto t_vvoo = btas::Tensor<double>{btas::Range{nvirt, nvirt, nocc, nocc}};
  t_vo.fill(0);
  t_vvoo.fill(0);

  auto d_vo = btas::Tensor<double>{btas::Range{nvirt, nocc}};
  auto d_vvoo = btas::Tensor<double>{btas::Range{nvirt, nvirt, nocc, nocc}};
  d_vo.fill(0.);
  d_vvoo.fill(0.);
  for (auto a = 0; a < nvirt; ++a)
    for (auto i = 0; i < nocc; ++i) {
      d_vo(a, i) = fock(i, i) - fock(nocc + a, nocc + a);
      for (auto b = 0; b < nvirt; ++b)
        for (auto j = 0; j < nocc; ++j)
          d_vvoo(a, b, i, j) =
              d_vo(a, i) + fock(j, j) - fock(nocc + b, nocc + b);
    }

  // ============= SeQuant ==================== //

  using sequant::eqs::cceqvec;
  using sequant::optimize::optimize;
  using sequant::optimize::tail_factor;
  using sequant::utils::binarize_expr;
  using evxpr_node = sequant::utils::binary_node<sequant::utils::eval_expr>;

  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::mbpt::set_default_convention();
  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  auto cc_r = cceqvec{2, 2}(true, true, true, true);

  auto r1_node = optimize(tail_factor(cc_r[1]));
  auto r2_node = optimize(tail_factor(cc_r[2]));

  auto eval_inst_r1 = eval_instance{r1_node};
  auto eval_inst_r2 = eval_instance{r2_node};

  auto yielder = yield_leaf{nocc, nvirt, fock, eri, t_vo, t_vvoo};

  auto manager = sequant::eval::cache_manager<btas::Tensor<double>>();

  auto visitor = [&manager](evxpr_node const& node) {
    if (node->tensor().label() == L"I")
      manager.add_key_decaying(node->hash());
    else if (node->tensor().label() != L"t")
      manager.add_key_persistent(node->hash());
    else
      return;  // tensor of type t_vo, t_vvoo
  };

  r1_node.visit(visitor);
  r2_node.visit(visitor);

  manager.keep_more_repeating();
  auto const decaying_imeds = manager.key_counts_decaying();

  auto const g_vvoo =
      yielder(sequant::utils::parse_expr(L"g_{a1,a2}^{i1,i2}",
                                         sequant::Symmetry::antisymm)
                  ->as<sequant::Tensor>());
  auto const f_vo = yielder(
      sequant::utils::parse_expr(L"f_{a1}^{i1}", sequant::Symmetry::antisymm)
          ->as<sequant::Tensor>());

  const auto maxiter = 100;
  const auto conv = 1e-12;

  size_t iter = 0;
  auto ediff = 0.0;
  auto normdiff = 0.0;
  auto ecc = 0.0;

  auto start = std::chrono::high_resolution_clock::now();
  auto empty_manager = sequant::eval::cache_manager<btas::Tensor<double>>();
  do {
    ++iter;
    auto r1 = eval_inst_r1.evaluate_asymm(yielder, manager);
    auto r2 = eval_inst_r2.evaluate_asymm(yielder, manager);

    auto norm_last = std::sqrt(btas::dot(t_vvoo, t_vvoo));

    // updating T1 and T2
    for (auto a = 0; a < nvirt; ++a)
      for (auto i = 0; i < nocc; ++i) {
        t_vo(a, i) += r1(a, i) / d_vo(a, i);
        for (auto b = 0; b < nvirt; ++b)
          for (auto j = 0; j < nocc; ++j)
            t_vvoo(a, b, i, j) += r2(a, b, i, j) / d_vvoo(a, b, i, j);
      }
    normdiff = norm_last - std::sqrt(btas::dot(t_vvoo, t_vvoo));

    // calculating ecc

    auto ecc_last = ecc;
    decltype(r2) temp;

    btas::contract(1.0, g_vvoo, {'a', 'b', 'i', 'j'}, t_vo, {'a', 'i'}, 0.0,
                   temp, {'b', 'j'});
    ecc = 0.5 * btas::dot(temp, t_vo)         //
          + 0.25 * btas::dot(g_vvoo, t_vvoo)  //
          + btas::dot(f_vo, t_vo);

    ediff = ecc_last - ecc;

    cout << "E(CC) = "
         << std::setprecision(std::numeric_limits<double>::max_digits10) << ecc
         << endl;

    manager.add_key_decaying(decaying_imeds);
  } while (iter < maxiter &&
           (std::fabs(normdiff) > conv || std::fabs(ediff) > conv));

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  cout << "\nOut of loop after " << iter << " iterations.\n"
       << "\nTime: " << duration.count() << " microseconds." << endl;

  return 0;
}

// testing symmetrization and anti-symmetrization
//
// auto node = binarize_expr(r1);
// auto eval_inst = evaluate_btas{node};
// auto eval_res = eval_inst.evaluate(yielder);
//
// cout << eval_res << " norm = " << norm(eval_res) << endl;

// t_vvoo.generate([]() { return static_cast<double>(std::rand()) / RAND_MAX;
// }); cout << "norm(t_vvoo) = " << norm(t_vvoo) << endl;
//
// auto expr = sequant::utils::parse_expr(L"t_{a1, a2}^{i1, i2}",
//                                        sequant::Symmetry::antisymm);
// auto node = binarize_expr(expr);
// auto eval_inst = evaluate_btas{node};
//
// // manual symmetrization
//
// auto manual_symm = decltype(t_vvoo){t_vvoo.range()};
// manual_symm.fill(0.);
//
// auto temp = manual_symm;
// btas::permute(t_vvoo, {0, 1, 2, 3}, temp, {0, 1, 2, 3});
// manual_symm += temp;
// btas::permute(t_vvoo, {0, 1, 2, 3}, temp, {1, 0, 3, 2});
// manual_symm += temp;
//
// cout << "norm(manual_symm) = " << norm(manual_symm) << endl;
// auto eval_symm = eval_inst.evaluate_symm(yielder);
// cout << "norm(eval_symm) = " << norm(eval_symm) << endl;
//
// // manual anti-symmetrization
// auto manual_asymm = decltype(t_vvoo){t_vvoo.range()};
// manual_asymm.fill(0.);
//
// temp = manual_asymm;
// btas::permute(t_vvoo, {0, 1, 2, 3}, temp, {0, 1, 2, 3});
// manual_asymm += temp;
// btas::permute(t_vvoo, {0, 1, 2, 3}, temp, {1, 0, 3, 2});
// manual_asymm += temp;
//
// btas::permute(t_vvoo, {0, 1, 2, 3}, temp, {1, 0, 2, 3});
// manual_asymm -= temp;
// btas::permute(t_vvoo, {0, 1, 2, 3}, temp, {0, 1, 3, 2});
// manual_asymm -= temp;
//
// cout << "norm(manual_asymm) = " << norm(manual_asymm) << endl;
// auto eval_asymm = eval_inst.evaluate_asymm(yielder);
// cout << "norm(eval_asymm) = " << norm(eval_asymm) << endl;
