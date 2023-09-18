#ifndef SEQUANT_ANALYSIS_UTILS_HPP
#define SEQUANT_ANALYSIS_UTILS_HPP

#include <SeQuant/core/expr.hpp>
#include <filesystem>

namespace sequant {

enum struct SpinTraceAlgorithm { Fast, Robust, Invalid, OpenShell = Invalid };

///
/// \brief Write line-by-line the result of \c sequant::deparse_expr on
///        each term of \c sum. std::exit(1) if the file cannot be written.
///
/// \param file The file to write to.
/// \param sum The sum of terms to write.
/// \warning This function will overwrite the file if it exists.
///
void write_sum(Sum const& sum, std::filesystem::path const& file);

void write_sum(std::ostream& os, Sum const& sum);

///
/// \brief Read a sum of terms from a file. std::exit(1) if file cannot be read.
///
/// \param file The file to read from.
/// \return ExprPtr to a sequant::Sum.
///
ExprPtr read_sum(std::filesystem::path const& file);

///
/// \brief Generate a canonical filename for coupled-cluster method expressed
///        in spinorbital or closedshell spinfree orbtials.
///
/// \param excit The highest level of excitation in the CC method.
/// \param r  The residual equation this filename will correspond to.
/// \param st_algo Spintrace algorithm type. Value equal to 'Invalid' implies
///                no spintrace.
/// \param optimize If true filename indicates equations are optimized.
/// \param ext Append this extension to the filename.
/// \return A canonical file name that can be useful to write/read to/from.
///
std::string canonical_fname_cck(size_t excit,                //
                                size_t r,                    //
                                SpinTraceAlgorithm st_algo,  //
                                bool optimize,               //
                                std::string_view ext);

///
/// \brief Generate a canonical filename for coupled-cluster method expressed
///        in openshell spinfree orbitals.
///
/// \param excit The highest level of excitation in the CC method.
/// \param r  The residual equation this filename will correspond to.
/// \param nalpha The number of particle indices with alpha tag.
///               0 <= \c nalpha <= r.
/// \param optimize If true filename indicates equations are optimized.
/// \param ext Append this extension to the filename.
/// \return A canonical file name that can be useful to write/read to/from.
///
std::string canonical_fname_cck(size_t excit,   //
                                size_t r,       //
                                size_t nalpha,  //
                                bool optimize,  //
                                std::string_view ext);

///
/// \brief Generate coupled-cluster equations in spinorbitals.
/// \param excit The highest level of excitation in the CC method.
/// \return A vector of ExprPtr to Sum.
/// \note The first tensor of each term in the sums (the OpType::A tensor) will
///       be removed in the returned expressions.
///
container::vector<ExprPtr> cck_equations(size_t excit);

///
/// \param excit The highest level of excitation in the CC method.
/// \param st_algo Spintracing algorithm to use. Fast, Robust or OpenShell.
/// \return A vector of ExprPtr to Sum.
/// \note The first tensor of each term in the sums(the OpType::S tensor) will
///       be removed in the returned expressions. This statement only applies
///       when \c st_algo is SpinTraceAlgorithm::{Fast,Robust}.
///
container::vector<ExprPtr> cck_equations(size_t excit,
                                         SpinTraceAlgorithm st_algo);

container::vector<std::string> canonical_fnames_cck(size_t excit, bool optimize,
                                                    std::string_view ext);

container::vector<std::string> canonical_fnames_cck(size_t excit,
                                                    SpinTraceAlgorithm st_algo,
                                                    bool optimize,
                                                    std::string_view ext);

inline static constexpr std::string_view database_dir_name = ".sequant";

std::filesystem::path database_dir();

void initialize_db();

}  // namespace sequant

#endif  // SEQUANT_ANALYSIS_UTILS_HPP
