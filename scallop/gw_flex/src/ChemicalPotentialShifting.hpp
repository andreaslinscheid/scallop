/*	This file ChemicalPotentialShifting.hpp is part of scallop.
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
 *  Created on: Dec 15, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/ChemicalPotentialShifting.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
void ChemicalPotentialShifting<T>::determine_new_chemPot(
		bT Nelectrons,
		bT invTemp,
		KohnShamBandStructure<T> const& elstr,
		KohnShamGreensFunctionOrbital<T> const& gks,
		GreensFunctionOrbital<T> const& g,
		std::pair<bool,bT> & newChemPot)
{
	assert( elstr.get_spaceGrid_proc().get_grid() == gks.get_spaceGrid_proc().get_grid() );
	assert( elstr.get_spaceGrid_proc().get_grid() == g.get_spaceGrid_proc().get_grid() );
	assert( gks.get_num_time() == g.get_num_time() );
	assert( gks.get_data_block_size() == g.get_data_block_size() );
	assert( gks.get_num_time() == g.get_num_time() );
	assert( (!gks.is_in_time_space()) && (!g.is_in_time_space()) );
	assert( (! gks.is_in_k_space()) && (! g.is_in_k_space()) );
	assert( elstr.get_nOrb() == g.get_nOrb() );

	size_t nM = g.get_num_time();
	size_t nB = g.get_data_block_size();
	size_t dataSize = nM*nB;
	if ( buffer2_.size() != dataSize)
		buffer2_ = buffer1_ = V(dataSize);

	MemoryLayout gflayout;
	gflayout.initialize_layout_2pt_obj( g.get_nOrb() );

	auxillary::LinearAlgebraInterface<T> linalg;

	auto compute_N_chem = [&] (bT deltaMu)
		{
			bT nElec = 0;
			assert( std::abs(elstr.get_chem_pot()-g.get_chem_pot())<0.00001 );
			//deltaMu is measured in absolute units, i.e. the band structure,
			// not relative to the current chemical potential.
			bT shift_mu = deltaMu - elstr.get_chem_pot();
			this->compute_N_elec_mB(
					gks,g,
					invTemp,shift_mu,nElec,
					gflayout);
			nElec += elstr.compute_N_electrons( invTemp, shift_mu );

			return (nElec-Nelectrons);
		};

	auto check_conv = [] (bT f1, bT f2)
		{
			return std::abs(f1-f2)<1e-4;
		};

	auto energyRange = elstr.get_band_width();
	energyRange.first += elstr.get_chem_pot();
	energyRange.second += elstr.get_chem_pot();

	auxillary::BasicFunctions bf;
	bf.secant_method( elstr.get_chem_pot(),
			energyRange,
			compute_N_chem,check_conv,1000,
			newChemPot);
}

template<typename T>
void ChemicalPotentialShifting<T>::compute_N_elec_mB(
		KohnShamGreensFunctionOrbital<T> const& gks,
		GreensFunctionOrbital<T> const& g,
		bT beta,
		bT muDelta,
		bT & nElec,
		MemoryLayout & gfl)
{
	g.chem_pot_adj_local_part( muDelta, buffer1_ );
	gks.chem_pot_adj_local_part( muDelta, buffer2_ );

	size_t nM = g.get_num_time();
	size_t nB = g.get_data_block_size();
	size_t nO = g.get_nOrb();

	T tmp = 0;
	for ( size_t iw = 0 ; iw < nM; ++iw )
		for ( size_t iO = 0; iO < nO; ++iO)
			for ( size_t a = 0; a < 2; ++a)
				for ( size_t s = 0; s < 2; ++s)
				tmp += (0.5/beta)*(buffer1_[iw*nB+gfl.memory_layout_2pt_obj(iO,a,s,iO,a,s)]
				           - buffer2_[iw*nB+gfl.memory_layout_2pt_obj(iO,a,s,iO,a,s)]);
	nElec += std::real(tmp);

	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	mpi.sum( nElec );
}

} /* namespace gw_flex */
} /* namespace scallop */
