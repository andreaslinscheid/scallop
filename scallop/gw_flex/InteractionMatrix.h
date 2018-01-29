/*	This file InteractionMatrix.h is part of scallop.
 *
 *  scallop is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  scallop is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with scallop.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Nov 22, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_INTERACTIONMATRIX_H_
#define SCALLOP_GW_FLEX_INTERACTIONMATRIX_H_

#include "scallop/auxillary/TemplateTypedefs.h"
#include "scallop/gw_flex/MemoryLayout.h"
#include <string>

namespace scallop
{
namespace gw_flex
{

template<typename T>
class InteractionMatrix : private MemoryLayout
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;

	using MemoryLayout::get_nOrb;
	using MemoryLayout::get_nChnls;

	void init_file( std::string const& filename );

	T * write_ptr();

	T const * read_ptr() const;

	T & operator() (size_t j, size_t jp, size_t l1, size_t l2, size_t l3, size_t l4);

	T operator() (size_t j, size_t jp, size_t l1, size_t l2, size_t l3, size_t l4) const;

	bool empty() const;
private:

	typename auxillary::TemplateTypedefs<T>::scallop_vector data_;

};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/InteractionMatrix.hpp"
#endif /* SCALLOP_GW_FLEX_INTERACTIONMATRIX_H_ */
