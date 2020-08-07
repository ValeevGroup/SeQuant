#include "SeQuant/core/attr.hpp"
#include "SeQuant/core/tensor_network.hpp"
#include "catch.hpp"

#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/factorize/inferred_imed_data.hpp>

using namespace sequant;
using namespace sequant::factorize;

TEST_CASE("Test InferredImedData", "[factorization][intermediate]") {
  SECTION("Product") {
    auto tnsr1 = std::make_shared<Tensor>(
        Tensor(L"g", {L"i_3", L"i_4"}, {L"i_1", L"i_2"}, Symmetry::antisymm));
    auto tnsr2 = std::make_shared<Tensor>(
        Tensor(L"t", {L"a_1", L"a_2"}, {L"i_3", L"i_4"}, Symmetry::antisymm));

    for (const auto& t : {tnsr1, tnsr2}) {
      assert(t->braket_symmetry() == BraKetSymmetry::conjugate);
      assert(t->symmetry() == Symmetry::antisymm);
      assert(t->particle_symmetry() == ParticleSymmetry::symm);
    }

    auto inferredData = InferredImedData(tnsr1, tnsr2);
    //
    // the whole bra of tnsr1 is contracted with the whole bra of tnsr2
    //
    REQUIRE_FALSE(inferredData.is_sum);
    REQUIRE(inferredData.symmetry == Symmetry::antisymm);
    REQUIRE(inferredData.particle_symmetry == ParticleSymmetry::symm);
    // REQUIRE(inferredData.braket_symmetry == BraKetSymmetry::conjugate);

    // contraction that doesn't allow to deduce full symmetry of the resulting
    // tensor
    auto tnsr3 = std::make_shared<Tensor>(
        Tensor(L"g", {L"i_3", L"a_1"}, {L"i_1", L"a_3"}, Symmetry::antisymm));
    auto tnsr4 = std::make_shared<Tensor>(
        Tensor(L"t", {L"a_2", L"a_3"}, {L"i_2", L"i_3"}, Symmetry::antisymm));
    inferredData = InferredImedData(tnsr3, tnsr4);
    REQUIRE_FALSE(inferredData.is_sum);
    REQUIRE(inferredData.symmetry == Symmetry::nonsymm);
    REQUIRE(inferredData.particle_symmetry == ParticleSymmetry::nonsymm);

    // outer product of the same tensor
    // Turned off for now.
    // auto tnsr5 = std::make_shared<Tensor>(Tensor(L"f", {L"i_1"}, {L"i_2"},
    // Symmetry::nonsymm)); auto tnsr6 = std::make_shared<Tensor>(Tensor(L"f",
    // {L"i_3"}, {L"i_4"}, Symmetry::nonsymm)); inferredData =
    // InferredImedData(tnsr5, tnsr6); REQUIRE_FALSE(inferredData.is_sum);
    // REQUIRE(inferredData.symmetry == Symmetry::symm);
    // REQUIRE(inferredData.particle_symmetry == ParticleSymmetry::symm);
  }

  SECTION("Sum") {
    auto tnsr1 = std::make_shared<Tensor>(
        Tensor(L"g", {L"a_1", L"a_2"}, {L"i_1", L"i_2"}, Symmetry::antisymm));
    auto tnsr2 = std::make_shared<Tensor>(
        Tensor(L"I", {L"a_1", L"a_2"}, {L"i_1", L"i_2"}, Symmetry::antisymm));

    auto inferredData = InferredImedData(tnsr1, tnsr2);

    REQUIRE(inferredData.is_sum);
    REQUIRE(inferredData.symmetry == Symmetry::antisymm);
    REQUIRE(inferredData.particle_symmetry == ParticleSymmetry::symm);
    // REQUIRE(inferredData.braket_symmetry == BraKetSymmetry::conjugate);

    auto tnsr3 = std::make_shared<Tensor>(
        Tensor(L"I", {L"a_1", L"a_2"}, {L"i_1", L"i_2"}, Symmetry::symm));
    inferredData = InferredImedData(tnsr1, tnsr3);
    REQUIRE(inferredData.is_sum);
    REQUIRE(inferredData.symmetry == Symmetry::symm);
    REQUIRE(inferredData.particle_symmetry == ParticleSymmetry::symm);

    auto tnsr4 = std::make_shared<Tensor>(
        Tensor(L"I", {L"a_1", L"a_2"}, {L"i_1", L"i_2"}, Symmetry::symm));

    inferredData = InferredImedData(tnsr3, tnsr4);
    REQUIRE(inferredData.is_sum);
    REQUIRE(inferredData.symmetry == Symmetry::symm);
    REQUIRE(inferredData.particle_symmetry == ParticleSymmetry::symm);
  }
}
