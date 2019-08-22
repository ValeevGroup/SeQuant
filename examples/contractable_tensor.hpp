#ifndef SEQUANT_CONTRACTABLE_TENSOR_HPP
#define SEQUANT_CONTRACTABLE_TENSOR_HPP

#include <btas/btas.h>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor.hpp>

using STensor   = sequant::Tensor;
using BTensor   = btas::Tensor<double>;
using index_vec = std::vector<sequant::Index>;

namespace sequant {

  namespace contractable {

    class Tensor
    {
      private:
        // using i_1, i_2, i_3, ... for occupied and
        //        a_1, a_2, a_3, ... for virtuals
        // using a sequant::Index object i_1 would be represented by Index{L"i_1"}
        //
        // NOTE: occ_index_ may have all or some virtual indices depending
        // on the type of the mathematical tensor it represents
        //
        // Characteristic to Tensor object is that if a virtual index appears
        // in occ_index_, then at the same position in virt_index_, another
        // virtual index must appear (an occupied index cannot appear in
        // virt_index_ if a virtual index appears at the equivalent position in occ_index_).
        // That is to say, the indices are sorted. The constructor takes care of
        // sorting the indices by keping the bra<->ket mapping.
        //
        index_vec occ_index_; // a vector of occupied indices
        index_vec virt_index_; // a vector of virtual indices
        // the tensor with data
        BTensor btensor_;
        // the label
        std::wstring label_ = L"";

      public:

        Tensor() = default;

        explicit Tensor(const STensor &st) {

          this->label_ = st.label();

          std::copy(st.bra().begin(), st.bra().end(), std::back_inserter(this->occ_index_));
          std::copy(st.ket().begin(), st.ket().end(), std::back_inserter(this->virt_index_));

          // sorting tensor Index objects
          // which should help us to map a given
          // tensor to the equivalent tensor that might have been already
          // computed
          const auto dim = st.rank();
          for (auto i = 0; i < dim; ++i) {
            if (st.bra()[i].label()[0] == st.ket()[i].label()[0])
              continue; // found oo or vv
            else {
              if (! (st.bra()[i] < st.ket()[i]) ) {
                std::swap(this->occ_index_, this->virt_index_);
              }
              break;
            }
          }
        }

        Tensor(const std::wstring& l, const index_vec &ov): Tensor() {
          //
          // if a vector of indices is passed to the constructor
          // the vector should be divisible into occupied half
          // and the virtual half spaces
          //
          auto size = ov.size();
          if (size%2 != 0)
            throw "The dimension of occupied and virtual space should match!";

          using std::back_inserter;
          using std::copy;
          // assumes the function parameter 'ov' is a sorted vector
          copy(ov.begin(), ov.begin()+size/2, back_inserter(this->occ_index_));
          copy(ov.begin()+size/2, ov.end(), back_inserter(this->virt_index_));
          //
          this->label_ = l;
        }

        Tensor operator+(const Tensor& other) {
          // the tensor with data pointed by btensor_ptr_ should have the same
          // shapes however, those conditions are NOT checked here
          Tensor result{other};
          result.btensor_ += this->btensor_;
          return result;
        }

        const index_vec&  occs() const { return this->occ_index_;  }

        const index_vec& virts() const { return this->virt_index_; }

        const std::wstring& label() const { return this->label_; }

        const BTensor& bt() const { return this->btensor_; }

        void link_btensor(const BTensor &bt) { this->btensor_ = bt; }

        // use with care
        void scale_btensor(const double& d) { btas::scal(d, this->btensor_); }

        std::wstring translate() const {
          std::wstring result = L"";
          auto ovify = [&result](const auto& idx_container){
            for (const auto& i: idx_container) {
              auto oclass = sequant::occupancy_class(i.space());
              if (oclass == -1)      // -1 = occupied
                result += 'o';
              else if (oclass == 1)  // 1 = virtuals
                result += 'v';
              else                   // 
                abort();
            }
          };
          ovify(occ_index_);
          ovify(virt_index_);
          return result;
        }
    };

  } // namespace sequant::contractable
} // namespace sequant

#endif /* SEQUANT_CONTRACTABLE_TENSOR_HPP */
