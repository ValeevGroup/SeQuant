#ifndef SEQUANT_EXTERNAL_INTERFACE_UTILS_HPP
#define SEQUANT_EXTERNAL_INTERFACE_UTILS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/result_expr.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <map>
#include <string>

class IndexSpaceMeta {
public:
	struct Entry {
		std::wstring tag;
		std::wstring name;
	};

	IndexSpaceMeta() = default;

	std::size_t getSize(const sequant::IndexSpace &space) const;

	std::size_t getSize(const sequant::Index &index) const;

	std::wstring getLabel(const sequant::IndexSpace &space) const;

	std::wstring getName(const sequant::IndexSpace &space) const;

	std::wstring getTag(const sequant::IndexSpace &space) const;

	void registerSpace(sequant::IndexSpace space, Entry entry);

private:
	std::map< sequant::IndexSpace, Entry > m_entries;
};

sequant::container::svector< sequant::container::svector< sequant::Index > >
	getExternalIndexPairs(const sequant::ExprPtr &expression);

bool needsSymmetrization(const sequant::ExprPtr &expression);

sequant::ExprPtr generateResultSymmetrization(const sequant::ResultExpr &result, std::wstring_view precursorName);

sequant::ExprPtr generateResultSymmetrization(const sequant::Tensor &result, std::wstring_view precursorName);

sequant::ExprPtr generateResultSymmetrization(std::wstring_view precursorName,
											  const sequant::IndexGroups< std::vector< sequant::Index > > &externals);

#endif
