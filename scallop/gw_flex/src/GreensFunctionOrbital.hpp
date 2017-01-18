/*	This file GreensFunctionOrbital.hpp is part of scallop.
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
 *  Created on: Nov 1, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/GreensFunctionOrbital.h"
#include "scallop/error_handling/error_handling.h"
#include "scallop/auxillary/LinearAlgebraInterface.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
GreensFunctionOrbital<T>::GreensFunctionOrbital()
	: MatsubaraImagTimeFourierTransform<T>(true)
{

}

template<typename T>
void GreensFunctionOrbital<T>::transform_itime_Mfreq_subtract(
		bT invT,
		GreensFunctionOrbital<T> const& ksGFTime,
		GreensFunctionOrbital<T> const& ksGFFreq )
{
	assert(  ksGFTime.is_in_time_space() );
	assert(  ! ksGFFreq.is_in_time_space() );
	assert( ksGFTime.get_spaceGrid_proc().get_grid() == this->get_spaceGrid_proc().get_grid() );
	assert( ksGFTime.get_num_time() == this->get_num_time() );
	assert( ksGFTime.get_data_block_size() == this->get_data_block_size() );
	assert( ksGFFreq.get_spaceGrid_proc().get_grid() == this->get_spaceGrid_proc().get_grid() );
	assert( ksGFFreq.get_num_time() == this->get_num_time() );
	assert( ksGFFreq.get_data_block_size() == this->get_data_block_size() );

	size_t nD = this->get_spaceGrid_proc().get_num_grid_data();
	size_t nM = this->get_num_time();
	size_t nB = this->get_data_block_size();

	T * ptrThisBegin = this->write_data_ptr_block(0,0);
	T const * ptrTimeBegin = ksGFTime.read_data_ptr_block(0,0);
	T const * ptrFreqBegin = ksGFFreq.read_data_ptr_block(0,0);

	if ( this->is_in_time_space() )
	{
		for ( size_t i = 0 ; i < nD*nM*nB ; ++i )
			ptrThisBegin[i] -= ptrTimeBegin[i];

		MatsubaraImagTimeFourierTransform<T>::transform_itime_Mfreq(invT);

		for ( size_t i = 0 ; i < nD*nM*nB ; ++i )
			ptrThisBegin[i] += ptrFreqBegin[i];
	}
	else // is in frequency space
	{
		for ( size_t i = 0 ; i < nD*nM*nB ; ++i )
			ptrThisBegin[i] -= ptrFreqBegin[i];

		MatsubaraImagTimeFourierTransform<T>::transform_itime_Mfreq(invT);

		for ( size_t i = 0 ; i < nD*nM*nB ; ++i )
			ptrThisBegin[i] += ptrTimeBegin[i];
	}
}

template<typename T>
void GreensFunctionOrbital<T>::initialize(
		size_t dimImTime,
		std::vector<size_t> gridDims,
		size_t orbitalDim,
		bool initialInTimeDomain,
		bool initialInReciprocalDomain,
		typename auxillary::TemplateTypedefs<T>::scallop_vector const& data,
		bT chemPot)
{
	chemPot_ = chemPot;

	this->initialize_layout_2pt_obj(orbitalDim);

	MatsubaraImagTimeFourierTransform<T>::initialize(
			dimImTime,
			gridDims,
			orbitalDim*orbitalDim*16,
			initialInTimeDomain,
			initialInReciprocalDomain,
			data);
}

template<typename T>
void GreensFunctionOrbital<T>::alter_chem_pot_no_shift( bT newChemicalPotential)
{
	chemPot_ = newChemicalPotential;
}

template<typename T>
void GreensFunctionOrbital<T>::set_chem_pot( bT newChemicalPotential)
{
	V alteredPart;
	assert( not this->is_in_k_space() );
	assert( not this->is_in_time_space() );
	this->chem_pot_adj_local_part(newChemicalPotential-chemPot_,alteredPart);
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	if ( mpi.get_mpi_me() == this->get_spaceGrid_proc().get_proc_index(/*in k space=*/false,/*iR=*/0) )
	{
		T * ptr = this->write_phs_grid_ptr_block(/*iR=*/0,/*iw*/0);
		std::copy( alteredPart.begin(), alteredPart.end(), ptr);
	}
	chemPot_ = newChemicalPotential;
}

template<typename T>
typename GreensFunctionOrbital<T>::bT
GreensFunctionOrbital<T>::get_chem_pot() const
{
	return chemPot_;
}

template<typename T>
T GreensFunctionOrbital<T>::operator() (
		size_t ig, size_t it, size_t l1, size_t a1, size_t s1,  size_t l2, size_t a2, size_t s2) const
{
	return *(this->read_phs_grid_ptr_block(ig,it) + this->memory_layout_2pt_obj(l1,a1,s1,l2,a2,s2) );
}

template<typename T>
T & GreensFunctionOrbital<T>::operator() (
		size_t ig, size_t it, size_t l1, size_t a1, size_t s1,  size_t l2, size_t a2, size_t s2)
{
	return *(this->write_phs_grid_ptr_block(ig,it) + this->memory_layout_2pt_obj(l1,a1,s1,l2,a2,s2) );
}

template<typename T>
T GreensFunctionOrbital<T>::operator() (
		size_t ig, size_t it, size_t m1,  size_t m2) const
{
	return *(this->read_phs_grid_ptr_block(ig,it) + this->memory_layout_combined_notation_2pt_obj(m1,m2) );
}

template<typename T>
T & GreensFunctionOrbital<T>::operator() (
		size_t ig, size_t it, size_t m1,  size_t m2)
{
	return *(this->write_phs_grid_ptr_block(ig,it) + this->memory_layout_combined_notation_2pt_obj(m1,m2) );
}

template<typename T>
void GreensFunctionOrbital<T>::chem_pot_adj_local_part(
		bT diffMu, V & localPart) const
{
	assert( (not this->is_in_k_space()) );

	size_t nB = this->get_data_block_size();
	if ( localPart.size() != nB * this->get_num_time())
		localPart = V( nB * this->get_num_time(), T(0) );

	//Only one processor has the local part.
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	if ( mpi.get_mpi_me() == this->get_spaceGrid_proc().get_proc_index( false, /*iR=*/0 ) )
	{

		size_t nO = this->get_nOrb();
		MemoryLayout gflayout;
		gflayout.initialize_layout_2pt_obj( nO );

		auxillary::LinearAlgebraInterface<T> linalg;

		bT dummy;
		V localBlock( nB );
		for ( size_t iw = 0 ; iw < this->get_num_time(); ++iw )
		{
			T const * ptr = this->read_phs_grid_ptr_block(/*iR=*/0,iw);
			std::copy( ptr, ptr+nB, localBlock.begin());
			linalg.invert_square_matrix(localBlock,dummy,false);
			for ( size_t iO = 0; iO < nO; ++iO)
				for ( size_t a = 0; a < 2; ++a)
					for ( size_t s = 0; s < 2; ++s)
						localBlock[gflayout.memory_layout_2pt_obj(iO,a,s,iO,a,s)] -=
							(a == 1? bT(-1.0):bT(1.0))*diffMu;
			linalg.invert_square_matrix(localBlock,dummy,false);
			std::copy( localBlock.begin(), localBlock.end(), localPart.begin()+iw*nB );
		}
	}
}

} /* namespace gw_flex */
} /* namespace scallop */
