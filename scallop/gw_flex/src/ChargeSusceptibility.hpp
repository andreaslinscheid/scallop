/*	This file ChargeSusceptibility.hpp is part of scallop.
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

#include "scallop/gw_flex/ChargeSusceptibility.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
void ChargeSusceptibility<T>::compute_from_gf( GreensFunctionOrbital<T> const & GF )
{
	GeneralizedSusceptibility<T>::compute_from_gf(GF,1);
}

template<typename T>
void ChargeSusceptibility<T>::copy_charge_part( GeneralizedSusceptibility<T> const& gsust )
{
	if ( not this->is_init() )
		this->initialize_zero( gsust.get_num_time(), gsust.get_nOrb(), 1,
				gsust.get_spaceGrid_proc().get_grid(), gsust.is_in_time_space(), gsust.is_in_k_space() );

	assert( this->get_nOrb() ==  gsust.get_nOrb() );
	assert( this->get_num_time() ==  gsust.get_num_time() );
	assert( this->get_spaceGrid_proc().get_grid() ==  gsust.get_spaceGrid_proc().get_grid() );
	assert( not (this->is_in_time_space() xor gsust.is_in_time_space()) );
	assert( not (this->is_in_k_space() xor gsust.is_in_k_space()) );

	size_t nG = this->is_in_k_space() ?
			this->get_spaceGrid_proc().get_num_k_grid() : this->get_spaceGrid_proc().get_num_R_grid();

	size_t dataBlockSize = std::pow(this->get_nOrb(),4);
	for ( size_t ig = 0; ig < nG; ++ig)
	{
		for ( size_t it = 0; it < this->get_num_time(); ++it)
		{
			T * blockPtr = this->write_phs_grid_ptr_block(ig,it);
			T const * blockPtrGSusct = gsust.read_phs_grid_ptr_block(ig,it);
			std::copy( blockPtrGSusct, blockPtrGSusct+dataBlockSize, blockPtr);
		}
	}
}

template<typename T>
void ChargeSusceptibility<T>::charge_RPA_enhancement(InteractionMatrix<T> const& interMat)
{
	this->spin_RPA_enhancement(interMat);
}

template<typename T>
T ChargeSusceptibility<T>::operator() (size_t ik, size_t iw, size_t m1,  size_t m2) const
{
	return GeneralizedSusceptibility<T>::operator ()(ik,iw,0,0,m1,m2);
}

template<typename T>
T & ChargeSusceptibility<T>::operator() (size_t ik, size_t iw, size_t m1,  size_t m2)
{
	return GeneralizedSusceptibility<T>::operator ()(ik,iw,0,0,m1,m2);
}

template<typename T>
T ChargeSusceptibility<T>::operator() (size_t ik, size_t iw, size_t l1, size_t l2, size_t l3, size_t l4) const
{
	return GeneralizedSusceptibility<T>::operator ()(ik,iw,0,0,l1,l2,l3,l4);
}

template<typename T>
T & ChargeSusceptibility<T>::operator() (size_t ik, size_t iw, size_t l1, size_t l2, size_t l3, size_t l4)
{
	return GeneralizedSusceptibility<T>::operator ()(ik,iw,0,0,l1,l2,l3,l4);
}

} /* namespace gw_flex */
} /* namespace scallop */
