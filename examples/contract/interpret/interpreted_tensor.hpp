//
// Created by Bimal Gaudel on 2019-09-08.
//

#ifndef SEQUANT_INTERPRETED_TENSOR_H
#define SEQUANT_INTERPRETED_TENSOR_H

#include <vector>
#include <string>
#include <algorithm>
#include <iterator>

#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor.hpp>
// TODO:  rewrite so that it requires only
// <SeQuant/core/expr_fwd.hpp>

namespace sequant {
  namespace interpret {
    template<typename T>
    class InterpretedTensor {
      public:
        InterpretedTensor<T>(const sequant::Tensor&);

        InterpretedTensor<T>(const std::wstring&, const std::vector<sequant::Index>&);

        const auto& label() const;

        const auto& bras() const;

        const auto& kets() const;

        const auto& tensor() const;

        auto translate() const;

        void link_tensor(const T&);

        // use with care
        void scale_tensor(const double&);

      protected:
        std::vector<sequant::Index> bra_indices_;

        std::vector<sequant::Index> ket_indices_;

        T tensor_;

        std::wstring label_ = L"";

      }; // class InterpretedTensor

    template<typename T>
      InterpretedTensor<T>::InterpretedTensor(const sequant::Tensor& st) {

        this->label_ = st.label();

        std::copy(st.bra().begin(), st.bra().end(),
            std::back_inserter(this->bra_indices_));

        std::copy(st.ket().begin(), st.ket().end(),
            std::back_inserter(this->ket_indices_));

        // sorting tensor Index objects
        // which should help us to map a given
        // tensor to the equivalent tensor that might have been already
        // computed
        //
        // NOTE: assumes antisymmetrized orbitals

        const auto dim = st.rank();
        for (auto i = 0; i < dim; ++i) {
          if (st.bra()[i].label()[0] == st.ket()[i].label()[0])
            continue; // found oo or vv
          else {
            if (! (st.bra()[i] < st.ket()[i]) ) {
              std::swap(this->bra_indices_, this->ket_indices_);
            }
            break;
          }
        }
      }

    template<typename T>
      InterpretedTensor<T>::InterpretedTensor(const std::wstring& lbl,
          const std::vector<sequant::Index>& ov) {

        // TODO: extend to handle different sized
        // bras and kets such as those encountered in EOM-CCSD
        auto size = ov.size();
        if (size%2 != 0)
          throw "The dimension of bras and kets should match!";

        using std::back_inserter;
        using std::copy;
        // assumes the function parameter 'ov' is a sorted vector
        copy(ov.begin(), ov.begin()+size/2,
            back_inserter(this->bra_indices_));
        copy(ov.begin()+size/2, ov.end(),
            back_inserter(this->ket_indices_));
        //
        this->label_ = lbl;
      }

    template <typename T>
      const auto& InterpretedTensor<T>::label() const {
        return this->label_;
      }

    template<typename T>
      const auto& InterpretedTensor<T>::bras() const {
        return this->bra_indices_;
      }

    template<typename T>
      const auto& InterpretedTensor<T>::kets() const {
        return this->ket_indices_;
      }

    template<typename T>
      const auto& InterpretedTensor<T>::tensor() const {
        return this->tensor_;
      }

    template<typename T>
      auto InterpretedTensor<T>::translate() const {

        std::wstring result = L"";
        auto ovify = [&result](const auto& idx_container){
          for (const auto& i: idx_container) {
            auto oclass = sequant::occupancy_class(i.space());
            if (oclass == -1)      // -1 = occupied
              result += 'o';
            else if (oclass == 1)  // 1 = virtuals
              result += 'v';
            else                   // 
              throw "Neither occupied nor unoccupied orbital encountered";
          }
        };

        ovify(this->bra_indices_);
        ovify(this->ket_indices_);
        return result;
      }

    template<typename T>
      void InterpretedTensor<T>::link_tensor(const T& t) {
        this->tensor_ = t;
      }

    // use with care
    template<typename T>
      void InterpretedTensor<T>::scale_tensor(const double& d) {
#ifdef SEQUANT_HAS_BTAS
        btas::scal(d, this->tensor_);
#endif
#ifdef SEQUANT_HAS_TILEDARRAY
        TA::scale(this->tensor_, d);
#endif
      }

  } // namespace interpret

}   // namespace sequant

#endif /* SEQUANT_INTERPRETED_TENSOR_H */
