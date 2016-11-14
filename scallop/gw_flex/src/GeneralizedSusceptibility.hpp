/*	This file GeneralizedSusceptibility.hpp is part of scallop.
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
 *  Created on: Nov 11, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/GeneralizedSusceptibility.h"
#include <assert.h>

namespace scallop
{
namespace gw_flex
{

template<typename T>
GeneralizedSusceptibility<T>::GeneralizedSusceptibility() :
	MatsubaraImagTimeFourierTransform<T>( /* bool Fermi=*/ false )
{

}

template<typename T>
T & GeneralizedSusceptibility<T>::operator() (size_t ik, size_t iw, size_t j, size_t jp, size_t m1,  size_t m2)
{
	return *(this->write_phs_grid_ptr_block(ik,iw)+this->memory_layout_combined_notation_4pt_scalar_obj(j,jp,m1,m2));
}

template<typename T>
T GeneralizedSusceptibility<T>::operator() (size_t ik, size_t iw, size_t j, size_t jp, size_t m1,  size_t m2) const
{
	return *(this->read_phs_grid_ptr_block(ik,iw)+this->memory_layout_combined_notation_4pt_scalar_obj(j,jp,m1,m2));
}

template<typename T>
void GeneralizedSusceptibility<T>::set_uninitialized()
{
	isInit_ = false;
}

template<typename T>
void GeneralizedSusceptibility<T>::compute_from_gf(GreensFunctionOrbital<T> const & GF)
{
	//check if we need to update, or initialize
	if ( ! isInit_ )
	{
		MemoryLayout::initialize_layout_4pt_scalar_obj( GF.get_nOrb() ,4);
		size_t nC = this->get_nChnls();

		size_t nM = GF.get_num_time();
		auto grid =  GF.get_spaceGrid_proc().get_grid();
		size_t blocksize = std::pow( GF.get_nOrb() , 4 )*nC*nC;
		typename auxillary::TemplateTypedefs<T>::scallop_vector data;

		// initialize an object with data being zero
		this->initialize( nM, grid, blocksize, /*In time space=*/ true, /*In k space=*/ true, data );
		isInit_ = true;
		bufferQ1_ = bufferQ2_ = typename auxillary::TemplateTypedefs<T>::scallop_vector( 4*GF.get_nOrb()*GF.get_nOrb()* nC*nC );
		bufferSBlock_ = typename auxillary::TemplateTypedefs<T>::scallop_vector( 16*std::pow(GF.get_nOrb(),4) );
	}

	size_t nO = this->get_nOrb();
	size_t nC = this->get_nChnls();

	assert( nO ==  GF.get_nOrb() );
	assert( this->get_num_time() ==  GF.get_num_time() );
	assert( this->get_spaceGrid_proc().get_grid() ==  GF.get_spaceGrid_proc().get_grid() );
	assert( (not GF.is_in_k_space()) and GF.is_in_time_space() );

	auxillary::LinearAlgebraInterface<T> linalg;

	MemoryLayout gf_layout;
	gf_layout.initialize_layout_2pt_obj( nO );
	/*
	 * NOTE: (TODO)
	 *
	 * We do not implement non-inversionally symmetric systems at this point.
	 * This would open a can of worms, since the inverse index may well be on another processor.
	 * Similarly, for non-local interactions, where this routine would have to be called with an R-offset.
	 * If both is not the case, then this is a GF local product (in space)!
	 */
	for ( size_t iR = 0; iR < this->get_spaceGrid_proc().get_num_R_grid(); ++iR)
	{
		for ( size_t it = 0; it < this->get_num_time(); ++it)
		{
			//Note that the shift here (-\tau) -> (beta - \tau) causes a minus sign to appear because
			// of the Fermionic Matsubara rules.
			size_t it_inverse = static_cast<size_t>(static_cast<int>( this->get_num_time() ) - static_cast<int>(it) - 1 );

			T * blockPtrSusct = this->write_phs_grid_ptr_block(iR,it);
			T const * blockPtrG1 = GF.read_phs_grid_ptr_block(iR,it_inverse);
			T const * blockPtrG2 = GF.read_phs_grid_ptr_block(iR,it);

			//Construct the 4Norbital x 16 matrix Q1 and the 16 x 4Norbital matrix Q2
			//We use a hard-coded matrix multiplication v_j * G
			for ( size_t j = 0 ; j < nC; ++j )
				for ( size_t l1 = 0 ; l1 < nO; ++l1 )
					for ( size_t a1 = 0 ; a1 < 2; ++a1 )
						for ( size_t s1 = 0 ; s1 < 2; ++s1 )
							for ( size_t l2 = 0 ; l2 < nO; ++l2 )
								for ( size_t a2 = 0 ; a2 < 2; ++a2 )
									for ( size_t s2 = 0 ; s2 < 2; ++s2 )
									{
										T prefactorG1;
										size_t gf1_index_NS_transpose;
										this->v_matrix_multiplication( j, l1, a2, s2, l2, a1, s1,
												gf_layout,gf1_index_NS_transpose ,prefactorG1);

										size_t Q1_index = (((((j*nO+l1)*nO+l2)*2+a1)*2+s1)*2+a2)*2+s2;
										assert( Q1_index < bufferQ1_.size() );
										bufferQ1_[Q1_index] = prefactorG1*blockPtrG1[ gf1_index_NS_transpose ];

										T prefactorG2;
										size_t gf2_index_NS;
										this->v_matrix_multiplication( j, l1, a1, s1, l2, a2, s2,
												gf_layout,gf2_index_NS ,prefactorG2);

										size_t Q2_index = (((((a1*2+s1)*2+a2)*2+s2)*4+j)*nO+l1)*nO+l2;
										assert( Q2_index < bufferQ2_.size() );
										bufferQ2_[Q2_index] =  prefactorG2*blockPtrG2[ gf2_index_NS ];
									}

			int dimSust = static_cast<int>(4*nO*nO);
			linalg.call_gemm(false, false,
					dimSust,dimSust,16,
					T(-1.0/16),	bufferQ1_.data(), 16,
								bufferQ2_.data(), dimSust,
		            T(0.0),bufferSBlock_.data(),dimSust);

			//Now, the susceptibility contains the right data, but not in the right order.
			//The current order is j,l1,l3,jp,l4,l2 if the first GF has orbital indices l1,l3 and the second l4,l2
			//We rearrange this to j,jp,l1,l2,l3,l4 in the previous nomenclature
			for ( size_t j = 0 ; j < nC; ++j )
				for ( size_t l1 = 0 ; l1 < nO; ++l1 )
					for ( size_t l3 = 0 ; l3 < nO; ++l3 )
						for ( size_t jp = 0 ; jp < nC; ++jp )
							for ( size_t l4 = 0 ; l4 < nO; ++l4 )
								for ( size_t l2 = 0 ; l2 < nO; ++l2 )
								{
									size_t bufferIndex = ((((j*nO+l1)*nO+l3)*4+jp)*nO+l4)*nO+l2;
									size_t susctIndex = ((((j*4+jp)*nO+l1)*nO+l2)*nO+l3)*nO+l4;
									blockPtrSusct[susctIndex] = bufferSBlock_[bufferIndex];
								}
		}
	}
}

template<typename T>
void GeneralizedSusceptibility<T>::v_matrix_multiplication(
		size_t j, size_t l1, size_t a1, size_t s1, size_t l2,  size_t a2, size_t s2,
		MemoryLayout const& gf_layout,
		size_t &gf_index, T & prefactor) const
{
	auto flip_spin = [] (size_t s) { return s == 0 ? 1 : 0; };

	switch (j) {
		case 0:
			prefactor = a1 != 0 ? T(1) : T(-1);
			gf_index = gf_layout.memory_layout_2pt_obj(l1,a1,s1,l2,a2,s2);
			break;
		case 1:
			prefactor = a1 != 0 ? T(1) : T(-1);
			gf_index = gf_layout.memory_layout_2pt_obj(l1,a1,flip_spin(s1),l2,a2,s2);
			break;
		case 2:
			prefactor = s1 != 0 ? T(0,1) : T(0,-1);
			gf_index = gf_layout.memory_layout_2pt_obj(l1,a1,flip_spin(s1),l2,a2,s2);
			break;
		case 3:
			prefactor = (s1 + a1 == 1) ? T(-1) : T(1);
			gf_index = gf_layout.memory_layout_2pt_obj(l1,a1,s1,l2,a2,s2);
			break;
		default:
			break;
	}
}

} /* namespace gw_flex */
} /* namespace scallop */
