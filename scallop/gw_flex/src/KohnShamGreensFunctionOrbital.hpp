/*	This file KohnShamGreensFunctionOrbital.hpp is part of scallop.
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

#include "scallop/gw_flex/KohnShamGreensFunctionOrbital.h"
#include "scallop/auxillary/LinearAlgebraInterface.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
void KohnShamGreensFunctionOrbital<T>::set_in_time_space(
		UnitaryWannierKSBands<T> unitaryWannierBands,
		std::vector<bT> KSBands,
		size_t timeDim,
		bT invTemp)
{
	unitary_ = std::move(unitaryWannierBands);
	KSBands_ = std::move(KSBands);
	this->set_in_both_spaces(timeDim,invTemp,true);
}

template<typename T>
void KohnShamGreensFunctionOrbital<T>::set_in_frequency_space(
		UnitaryWannierKSBands<T> unitaryWannierBands,
		std::vector<bT> KSBands,
		size_t freqDim,
		bT invTemp)
{
	unitary_ = std::move(unitaryWannierBands);
	KSBands_ = std::move(KSBands);
	this->set_in_both_spaces(freqDim,invTemp,false);
}

template<typename T>
void KohnShamGreensFunctionOrbital<T>::set_in_both_spaces(
		size_t timeDim,
		bT invTemp,
		bool timeSpace)
{

	const size_t nK = unitary_.get_spaceGrid_proc().get_num_k_grid();
	const size_t nM = timeDim;
	const size_t nO = unitary_.get_nOrb();

	auto fullGrid = unitary_.get_spaceGrid_proc().get_grid();

	if ( ! this->is_init() )
	{
		typename auxillary::TemplateTypedefs<T>::scallop_vector data;
		this->initialize(
				timeDim, fullGrid , nO, timeSpace, true, data);
	}

	auto FermiFunc = [] (bT E, bT beta)
		{
			return 1.0/ ( std::exp(E*beta) + 1.0 );
		};

	auxillary::LinearAlgebraInterface<T> linalg;
	std::vector<T> bareGF(nO*4);
	std::vector<T> conjUnitary(nO*4*nO*4);
	std::vector<T> tmp(nO*4*nO*4);
	for ( size_t ik = 0 ; ik < nK ; ++ik)
	{
		for ( size_t m1 = 0 ; m1 < nO*4 ; ++m1)
			for ( size_t m2 = 0 ; m2 < nO*4 ; ++m2)
				conjUnitary[m1*nO*4+m2] = std::conj(unitary_(ik,m2,m1));

		for ( size_t it = 0 ; it < nM ; ++it)
		{
			//Produce (-f(-E) exp(-E tau ))*A as a matrix in band/Nambu/spin space
			// since the above brackets is a diagonal matrix, we simply scale each row of A
			// by this matrix element.
			//We put the result already into this object and multiply the other unitary transform
			if( timeSpace )
			{
				for ( size_t iMKSB = 0 ; iMKSB < nO*4 ; ++iMKSB)
				{
					bT e = KSBands_[ik*nO*4+iMKSB];
					bT taui = invTemp*(bT(it)/bT(nM));
					bareGF[iMKSB] = -FermiFunc(-e,invTemp)*std::exp(-taui*e);
				}
			}
			else
			{
				for ( size_t iMKSB = 0 ; iMKSB < nO*4 ; ++iMKSB)
				{
					bT e = KSBands_[ik*nO*4+iMKSB];
					int frequencyIndex =
							(it < nM/2 ? static_cast<int>(it) : static_cast<int>(it)-static_cast<int>(nM) );

					bareGF[iMKSB] = 1.0 / (T(0,M_PI*(2*frequencyIndex+1)/invTemp) - e);
				}
			}
			auto unitaryThisK = unitary_.read_phs_grid_ptr_block(ik);
			linalg.matrix_times_diagonal_matrix( unitaryThisK, nO*4, bareGF.data(), tmp.data() );

			auto iteratorThis = this->write_phs_grid_ptr_block(ik,it);
			linalg.matrix_times_matrix( conjUnitary.data(), nO*4, tmp.data(), iteratorThis );
		}
	}
}

} /* namespace gw_flex */
} /* namespace scallop */
