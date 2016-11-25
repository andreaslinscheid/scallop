/*	This file SelfEnergy.hpp is part of scallop.
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
 *  Created on: Nov 24, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/SelfEnergy.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
void SelfEnergy<T>::set_to_zero()
{
	size_t dataDim = this->get_spaceGrid_proc().get_num_grid_data()
					*this->get_num_time()
					*std::pow( 4*this->get_nOrb()*this->get_nOrb(),2);
	std::fill( this->write_data_ptr_block(0,0),this->write_data_ptr_block(0,0)+dataDim, T(0) );
}

template<typename T>
void SelfEnergy<T>::add_electronic_selfenergy(
		GreensFunctionOrbital<T> const& gf,
		GeneralizedSusceptibility<T> const& sf)
{
	assert( gf.is_in_time_space() && sf.is_in_time_space() );
	assert( gf.get_num_time() == sf.get_num_time() );
	assert( (!gf.is_in_k_space()) && (!sf.is_in_k_space()) );
	assert( gf.get_nOrb() == sf.get_nOrb() );
	assert( gf.get_spaceGrid_proc().get_grid() == sf.get_spaceGrid_proc().get_grid() );

	if ( not this->is_init() )
	{
		this->initialize( gf.get_num_time(), gf.get_spaceGrid_proc().get_grid(), gf.get_nOrb(), true, false,
				typename auxillary::TemplateTypedefs<T>::scallop_vector() );
	}

	size_t nO = this->get_nOrb();
	size_t nM = this->get_num_time();
	size_t nC = sf.get_nChnls();

	assert( this->is_in_time_space() );
	assert( nM == sf.get_num_time() );
	assert( !this->is_in_k_space() );
	assert( nO == sf.get_nOrb() );
	assert( this->get_spaceGrid_proc().get_grid() == sf.get_spaceGrid_proc().get_grid() );

	typename auxillary::TemplateTypedefs<T>::scallop_vector bufferGF( 4*nO*4*nO );

	MemoryLayout gf_layout, sf_layout;
	gf_layout.initialize_layout_2pt_obj( nO );
	sf_layout.initialize_layout_4pt_scalar_obj( nO, nC );

	for ( size_t iR = 0; iR < sf.get_spaceGrid_proc().get_num_R_grid(); ++iR)
	{
		for ( size_t it = 0; it < sf.get_num_time(); ++it)
		{
			T * blockPtrSE = this->write_phs_grid_ptr_block(iR,it);
			T const * blockPtrGF = gf.read_phs_grid_ptr_block(iR,it);
			T const * blockPtrSust = sf.read_phs_grid_ptr_block(iR,it);

			for ( size_t j = 0 ; j < nC; ++j )
				for ( size_t jp = 0 ; jp < nC; ++jp )
				{
					this->transform_2pto_v(j,jp,nO,gf_layout,blockPtrGF,bufferGF.data());

					for ( size_t a1 = 0 ; a1 < 2; ++a1 )
						for ( size_t s1 = 0 ; s1 < 2; ++s1 )
							for ( size_t a2 = 0 ; a2 < 2; ++a2 )
								for ( size_t s2 = 0 ; s2 < 2; ++s2 )
									for ( size_t l1 = 0 ; l1 < nO; ++l1 )
										for ( size_t l2 = 0 ; l2 < nO; ++l2 )
											for ( size_t l3 = 0 ; l3 < nO; ++l3 )
												for ( size_t l4 = 0 ; l4 < nO; ++l4 )
												{
													size_t SEIndex = gf_layout.memory_layout_2pt_obj(l1,a1,s1,l2,a2,s2);
													size_t GFIndex = gf_layout.memory_layout_2pt_obj(l3,a1,s1,l4,a2,s2);
													size_t SustIndex = sf_layout.memory_layout_4pt_scalar_obj(j,jp,l1,l3,l2,l4);
													blockPtrSE[SEIndex] += blockPtrSust[SustIndex]*bufferGF[GFIndex];
												}
				}
		}
	}
}

template<typename T>
void SelfEnergy<T>::transform_2pto_v(
		size_t j, size_t jp, size_t nOrb, MemoryLayout const& twoPtLayout,
		T const * old2PtObj,
		T * transformed2PtObj) const
{
	auto flip_spin = [] (size_t s) { return s == 0 ? 1 : 0; };
	auto multiply_v_from_left = [&] ( size_t j, size_t a1, size_t s1, size_t &s1t, T & prefactor)
	{
		switch (j)
		{
			case 0:
				prefactor = a1 != 0 ? T(1) : T(-1);
				s1t = s1;
				break;
			case 1:
				prefactor = a1 != 0 ? T(1) : T(-1);
				s1t = flip_spin(s1);
				break;
			case 2:
				prefactor = s1 != 0 ? T(0,1) : T(0,-1);
				s1t = flip_spin(s1);
				break;
			case 3:
				prefactor = (s1 + a1 == 1) ? T(-1) : T(1);
				s1t = s1;
				break;
			default:
				break;
		}
	};
	auto multiply_v_from_right = [&] ( size_t j, size_t a2, size_t s2, size_t &s2t, T & prefactor)
	{
		switch (j)
		{
			case 0:
				prefactor = a2 != 0 ? T(1) : T(-1);
				s2t = s2;
				break;
			case 1:
				prefactor = a2 != 0 ? T(1) : T(-1);
				s2t = flip_spin(s2);
				break;
			case 2:
				prefactor = s2 != 1 ? T(0,1) : T(0,-1);
				s2t = flip_spin(s2);
				break;
			case 3:
				prefactor = (s2 + a2 == 1) ? T(-1) : T(1);
				s2t = s2;
				break;
			default:
				break;
		}
	};

	size_t s1t, s2t;
	T prefactor1, prefactor2;
	for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
		for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
		{
			multiply_v_from_left(j,a1,s1,s1t,prefactor1);
			for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
				for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
				{
					multiply_v_from_right(jp,a2,s2,s2t,prefactor2);
					for ( size_t l1 = 0 ; l1 < nOrb ; ++l1 )
						for ( size_t l2 = 0 ; l2 < nOrb ; ++l2 )
						{
							size_t iOld = twoPtLayout.memory_layout_2pt_obj(l1,a1,s1,l2,a2,s2);
							size_t iNew = twoPtLayout.memory_layout_2pt_obj(l1,a1,s1t,l2,a2,s2t);

							transformed2PtObj[iNew] = prefactor1*prefactor2*old2PtObj[iOld];
						}

				}
		}
};

} /* namespace gw_flex */
} /* namespace scallop */
