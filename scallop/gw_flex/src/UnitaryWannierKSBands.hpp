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
void UnitaryWannierKSBands<T>::initialize_identity(
		std::vector<size_t> spaceGrid,
		size_t numOrbitals)
{
	this->initialize_layout_2pt_obj(numOrbitals);
	size_t nO = this->get_nOrb();

	gdistr_.distribute_grid( std::move(spaceGrid) );
	size_t nK = this->get_spaceGrid_proc().get_num_k_grid();
	data_ = typename auxillary::TemplateTypedefs<T>::scallop_vector(16*nK*nO*nO);
	for ( size_t ik = 0 ; ik < nK ; ++ik)
		for ( size_t m1 = 0 ; m1 < nO*4 ; ++m1)
		{
			(*this)(ik,m1,m1) = T(1.0);
		}
}

template<typename T>
void UnitaryWannierKSBands<T>::initialize(
		std::vector<size_t> spaceGrid,
		size_t numOrbitals,
		typename auxillary::TemplateTypedefs<T>::scallop_vector data)
{
	this->initialize_layout_2pt_obj(numOrbitals);

	gdistr_.distribute_grid( std::move(spaceGrid) );
	data_ = std::move( data );
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
	return this->read_phs_grid_ptr_block(ik)+this->memory_layout_combined_notation_2pt_obj(m1,m2);
}

template<typename T>
T * UnitaryWannierKSBands<T>::write_phs_grid_ptr(size_t ik, size_t m1, size_t m2 )
{
	return this->write_phs_grid_ptr_block(ik)+this->memory_layout_combined_notation_2pt_obj(m1,m2);
}

template<typename T>
T const * UnitaryWannierKSBands<T>::read_phs_grid_ptr_block(size_t ik) const
{
	assert(ik*this->get_nOrb()*this->get_nOrb() < data_.size() );
	return &data_[ik*16*this->get_nOrb()*this->get_nOrb()];
}

template<typename T>
T * UnitaryWannierKSBands<T>::write_phs_grid_ptr_block(size_t ik)
{
	assert(ik*this->get_nOrb()*this->get_nOrb() < data_.size() );
	return  &data_[ik*16*this->get_nOrb()*this->get_nOrb()];
}

template<typename T>
parallel::GridDistribution<T> const&
UnitaryWannierKSBands<T>::get_spaceGrid_proc() const
{
	return gdistr_;
}

} /* namespace gw_flex */
} /* namespace scallop */
