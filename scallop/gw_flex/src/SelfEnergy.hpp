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
	if ( not this->is_init() )
		return;
	size_t dataDim = this->get_spaceGrid_proc().get_num_grid_data()
					*this->get_num_time() *std::pow( 4*this->get_nOrb(),2);
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
	size_t nC = sf.get_nChnls();

	assert( this->is_in_time_space() );
	assert( this->get_num_time() == sf.get_num_time() );
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
													blockPtrSE[SEIndex] += (nC == 4 ? -1.0:1.0)*blockPtrSust[SustIndex]*bufferGF[GFIndex];
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

	size_t s1t = 0, s2t = 0;
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

template<typename T>
void SelfEnergy<T>::linear_interpolate_frequency( V const & previousGrid, V const & presentGrid )
{
	assert( not this->is_in_time_space() );
	size_t nB = this->get_data_block_size();
	size_t nG = this->is_in_k_space() ? this->get_spaceGrid_proc().get_num_k_grid() : this->get_spaceGrid_proc().get_num_R_grid();

	assert( previousGrid.size() > 0 );
	VbT imagPreviousGrid( previousGrid.size() );
	for ( size_t i = 0 ; i < previousGrid.size(); ++i)
	{
		imagPreviousGrid[i] = i < previousGrid.size()/2 ?
				std::imag(previousGrid[i+previousGrid.size()/2]) :
				std::imag(previousGrid[static_cast<int>(i)-static_cast<int>(previousGrid.size()/2)]);
	}
	assert( presentGrid.size() > 0 );
	VbT imagPresentGrid( presentGrid.size() );
	for ( size_t i = 0 ; i < presentGrid.size(); ++i)
	{
		imagPresentGrid[i] = i < presentGrid.size()/2 ?
				std::imag(presentGrid[i+presentGrid.size()/2]) :
				std::imag(presentGrid[static_cast<int>(i)-static_cast<int>(presentGrid.size()/2)]);
	}

	auto map_fourier_order_back = [] (size_t i, size_t dim) {
		size_t r = i < dim/2 ?
				static_cast<int>(i)+dim/2
			    :static_cast<int>(i)-static_cast<int>(dim/2);
		assert( r < dim );
		return r;
	};

	V newData( imagPresentGrid.size()*nB*nG );
	for ( size_t ig = 0 ; ig < nG; ++ig)
	{
		for ( size_t it = 0 ; it < imagPresentGrid.size(); ++it)
		{
			size_t index = map_fourier_order_back(it, imagPresentGrid.size() );
			auto newDataPtr = &( newData[(ig*imagPresentGrid.size()+index)*nB] );

			bT Mfreq = imagPresentGrid[ it ];
			auto itl = std::lower_bound(imagPreviousGrid.begin(),imagPreviousGrid.end(),Mfreq);
			auto itu = std::upper_bound(imagPreviousGrid.begin(),imagPreviousGrid.end(),Mfreq);

			//If Mfreq is smaller than the range covered by the previous grid range
			//we simply copy the data from the smallest point.
			//Also if Mfreq hit the begin, exactly, we also simply copy this element
			if ( itl == imagPreviousGrid.begin() )
			{
				auto oldDataPtr = this->read_phs_grid_ptr_block(ig,this->get_num_time()/2);
				std::copy(oldDataPtr,oldDataPtr+nB,newDataPtr);
				continue;
			}

			//If Mfreq is larger than the range covered by the previous grid range
			//we simply copy the data from the largest point.
			//Also if freq hit the end, exactly, we also simply copy this element
			if ( itu == imagPreviousGrid.end() )
			{
				auto oldDataPtr = this->read_phs_grid_ptr_block(ig, this->get_num_time()/2-1);
				std::copy(oldDataPtr,oldDataPtr+nB,newDataPtr);
				continue;
			}

			//If the Mfreq is actually part of the previous set and in the current range,
			//	lower and upper bound don't give the same element.
			if ( itl != itu )
			{
				size_t indexOldGrid = std::distance(imagPreviousGrid.begin(),itl);
				auto oldDataPtr = this->read_phs_grid_ptr_block(ig,indexOldGrid);
				std::copy(oldDataPtr,oldDataPtr+nB,newDataPtr);
				continue;
			}

			--itl;//OK since itl cannot be begin() at this point.
			size_t indexLower = std::distance(imagPreviousGrid.begin(),itl);
			size_t indexUpper = std::distance(imagPreviousGrid.begin(),itu);
			assert( indexLower < imagPreviousGrid.size() );
			assert( indexUpper < imagPreviousGrid.size() );
			bT x = (Mfreq - imagPreviousGrid[indexLower])/(imagPreviousGrid[indexUpper]-imagPreviousGrid[indexLower]);
			assert( (x >= 0) && (x < 1.0) );

			auto ptrDatal = this->read_phs_grid_ptr_block(ig,map_fourier_order_back(indexLower,imagPreviousGrid.size()));
			auto ptrDatau = this->read_phs_grid_ptr_block(ig,map_fourier_order_back(indexUpper,imagPreviousGrid.size()));
			for ( size_t i = 0 ; i < nB; ++i)
				newDataPtr[i] = ptrDatal[i] * (bT(1)-x) + ptrDatau[i] * x;
		}
	}

	MatsubaraImagTimeFourierTransform<T>::initialize( imagPresentGrid.size(),
			this->get_spaceGrid_proc().get_grid(),
			this->get_nOrb()*16,
			this->is_in_time_space(),
			this->is_in_k_space(),
			newData);
}

template<typename T>
template<class G>
void SelfEnergy<T>::break_u1_symmetry(G const& gap,
		size_t nO,
		std::vector<size_t> const& fullGrid,
		V const& matsubaraGrid)
{
	if ( not this->is_init() )
	{
		this->initialize( matsubaraGrid.size(),
				fullGrid,
				nO,
				/*in time space*/false,
				/*in k space */true,
				typename auxillary::TemplateTypedefs<T>::scallop_vector() );
	}

	assert( this->is_in_k_space() );
	assert( not this->is_in_time_space() );

	size_t const nK = this->get_spaceGrid_proc().get_num_k_grid();
	size_t const nM = this->get_num_time();

	MemoryLayout ml;
	ml.initialize_layout_2pt_obj(nO);

	auto kgridProc = this->get_spaceGrid_proc().get_grid();
	for (size_t ik = 0 ; ik < nK ; ++ik)
	{
		auto kIndices = this->get_spaceGrid_proc().k_conseq_local_to_xyz_total(ik);
		std::vector<double> k( kIndices.size() );
		for ( size_t ikx = 0 ; ikx < k.size(); ++ikx )
			k[ikx] = double(kIndices[ikx])/double(kgridProc[ikx]);
		for (size_t iw = 0 ; iw < nM ; ++iw)
		{
			auto ptr = this->write_phs_grid_ptr_block(ik,iw);
			double omegan = double(std::imag(matsubaraGrid[iw]));
			for ( size_t l1 = 0 ; l1 < nO; ++l1)
			{
				//set the Nambu off-diagonal part with the singlet channel
				for ( size_t a1 = 0 ; a1 < 2; ++a1)
					for ( size_t s1 = 0 ; s1 < 2; ++s1)
						for ( size_t a2 = 0 ; a2 < 2; ++a2)
							for ( size_t s2 = 0 ; s2 < 2; ++s2)
								ptr[ml.memory_layout_2pt_obj(l1,a1,s1,l1,a2,s2)]
								    = gap(k,omegan,0,l1)*auxillary::BasicFunctions::singlet_Re_channel(a1,s1,a2,s2);

				//add the Nambu off-diagonal part with the tx channel
				for ( size_t a1 = 0 ; a1 < 2; ++a1)
					for ( size_t s1 = 0 ; s1 < 2; ++s1)
						for ( size_t a2 = 0 ; a2 < 2; ++a2)
							for ( size_t s2 = 0 ; s2 < 2; ++s2)
								ptr[ml.memory_layout_2pt_obj(l1,a1,s1,l1,a2,s2)]
								    += gap(k,omegan,1,l1)*auxillary::BasicFunctions::triplet_x_Re_channel(a1,s1,a2,s2);

				//add the Nambu off-diagonal part with the ty channel
				for ( size_t a1 = 0 ; a1 < 2; ++a1)
					for ( size_t s1 = 0 ; s1 < 2; ++s1)
						for ( size_t a2 = 0 ; a2 < 2; ++a2)
							for ( size_t s2 = 0 ; s2 < 2; ++s2)
								ptr[ml.memory_layout_2pt_obj(l1,a1,s1,l1,a2,s2)]
								    += gap(k,omegan,2,l1)*auxillary::BasicFunctions::triplet_y_Re_channel(a1,s1,a2,s2);

				//set the Nambu off-diagonal part with the tz channel
				for ( size_t a1 = 0 ; a1 < 2; ++a1)
					for ( size_t s1 = 0 ; s1 < 2; ++s1)
						for ( size_t a2 = 0 ; a2 < 2; ++a2)
							for ( size_t s2 = 0 ; s2 < 2; ++s2)
								ptr[ml.memory_layout_2pt_obj(l1,a1,s1,l1,a2,s2)]
								    += gap(k,omegan,3,l1)*auxillary::BasicFunctions::triplet_z_Re_channel(a1,s1,a2,s2);
			}
		}
	}
}

} /* namespace gw_flex */
} /* namespace scallop */
