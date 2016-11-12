/*	This file UnitaryWannierKSBands.hpp is part of scallop.
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
 *  Created on: Nov 3, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/UnitaryWannierKSBands.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
UnitaryWannierKSBands<T>::UnitaryWannierKSBands()
{

}

template<typename T>
void UnitaryWannierKSBands<T>::initialize_identity(
		std::vector<size_t> spaceGrid,
		size_t numOrbitals)
{
	this->initialize_layout_2pt_obj(numOrbitals);
	size_t nO = this->get_nOrb();

	typename auxillary::TemplateTypedefs<T>::scallop_vector data;
	this->initialize( data, true, false, std::move(spaceGrid), 1, 16*nO*nO );

	//Since we are initializing the identity everywhere, the grid does not matter.
	size_t nG = this->get_spaceGrid_proc().get_num_k_grid();

	for ( size_t ik = 0 ; ik < nG ; ++ik)
		for ( size_t m1 = 0 ; m1 < nO*4 ; ++m1)
		{
			(*this)(ik,m1,m1) = T(1.0);
		}
}

template<typename T>
T UnitaryWannierKSBands<T>::operator() (size_t ik, size_t m1, size_t m2) const
{
	return *(this->read_phs_grid_ptr(ik,m1,m2));
}

template<typename T>
T & UnitaryWannierKSBands<T>::operator() (size_t ik, size_t m1, size_t m2)
{
	return *(this->write_phs_grid_ptr(ik,m1,m2));
}

template<typename T>
T const * UnitaryWannierKSBands<T>::read_phs_grid_ptr(size_t ik, size_t m1, size_t m2 ) const
{
	return FFTBase<T>::read_phs_grid_ptr_block(ik,0)+this->memory_layout_combined_notation_2pt_obj(m1,m2);
}

template<typename T>
T * UnitaryWannierKSBands<T>::write_phs_grid_ptr(size_t ik, size_t m1, size_t m2 )
{
	return FFTBase<T>::write_phs_grid_ptr_block(ik,0)+this->memory_layout_combined_notation_2pt_obj(m1,m2);
}

template<typename T>
T const * UnitaryWannierKSBands<T>::read_phs_grid_ptr_block(size_t ik) const
{
	return FFTBase<T>::read_phs_grid_ptr_block(ik,0);
}

template<typename T>
T * UnitaryWannierKSBands<T>::write_phs_grid_ptr_block(size_t ik)
{
	return FFTBase<T>::write_phs_grid_ptr_block(ik,0);
}

} /* namespace gw_flex */
} /* namespace scallop */
