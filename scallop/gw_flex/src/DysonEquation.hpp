/*	This file DysonEquation.hpp is part of scallop.
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
 *  Created on: Dec 16, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/DysonEquation.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
void DysonEquation::solve_by_inversion(
		GreensFunctionOrbital<T> & g,
		typename auxillary::TemplateTypedefs<T>::scallop_vector const& MatsFreqs,
		KohnShamBandStructure<T> const& elstr,
		SelfEnergy<T> const& sigma) const
{
	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector V;
	typedef typename auxillary::TypeMapComplex<T>::type bT;

	size_t nK = sigma.get_spaceGrid_proc().get_num_k_grid();
	size_t nM = MatsFreqs.size();
	size_t nO = g.get_nOrb();
	size_t nOns = g.get_nOrb()*4;
	size_t nB = g.get_data_block_size();

	assert( sigma.get_spaceGrid_proc().get_grid() == g.get_spaceGrid_proc().get_grid());
	assert( elstr.get_spaceGrid_proc().get_grid() == sigma.get_spaceGrid_proc().get_grid());
	assert( nM == sigma.get_num_time());
	assert( nM == g.get_num_time());
	assert( sigma.is_in_k_space() && (!sigma.is_in_time_space()) );
	assert( sigma.get_data_block_size() == nB );
	assert( sigma.get_nOrb() == nO );
	assert( elstr.get_nOrb() == nO );

	g.set_time_space( false );
	g.set_k_space( true );

	MemoryLayout m;
	m.initialize_layout_2pt_obj( nO );

	auto w1 = [] (size_t spinNambu1, size_t spinNambu2)
		{
			return (spinNambu1 == spinNambu2 ? T(1.0) : T(0.0));
		};

	V KSHamOrb(nB);
	V tmp(nB);
	V conjUnitary(nB);
	V bnd(nOns);
	bT dummy = bT(0);
	auxillary::LinearAlgebraInterface<T> linalg;
	for ( size_t ik = 0 ; ik < nK; ++ik)
	{
		//Construct the KS band structure in orbital space.
		auto unitaryk = elstr.get_unitary().read_phs_grid_ptr_block(ik);

		for ( size_t l1 = 0 ;l1 < nO; ++l1)
			for ( size_t ns1 = 0 ;ns1 < 4; ++ns1)
				for ( size_t l2 = 0 ;l2 < nO; ++l2)
					for ( size_t ns2 = 0 ;ns2 < 4; ++ns2)
						conjUnitary[m.memory_layout_2pt_obj_nsc(l1,ns1,l2,ns2)]
						            = std::conj(unitaryk[m.memory_layout_2pt_obj_nsc(l2,ns2,l1,ns1)]);

		for ( size_t n = 0 ; n < nOns; ++n)
			bnd[n] = elstr.get_bands()[ik*nOns+n];
		linalg.matrix_times_diagonal_matrix( unitaryk, nOns, bnd.data(), tmp.data() );
		linalg.matrix_times_matrix( conjUnitary.data(), nOns, tmp.data(), KSHamOrb.data() );

		//This constructs the inverse Green's function and inverts the matrix in orbital space
		// G = (1 - H^KS - Simga)^-1
		for ( size_t iw = 0 ; iw < nM; ++iw)
		{
			auto sePtr = sigma.read_phs_grid_ptr_block(ik,iw);
			for ( size_t l1 = 0 ;l1 < nO; ++l1)
				for ( size_t ns1 = 0 ;ns1 < 4; ++ns1)
					for ( size_t l2 = 0 ;l2 < nO; ++l2)
						for ( size_t ns2 = 0 ;ns2 < 4; ++ns2)
						{
							size_t index =  m.memory_layout_2pt_obj_nsc(l1,ns1,l2,ns2);
							tmp[index] = bT(l1==l2?1.0:0.0)*w1(ns1,ns2)*MatsFreqs[iw] - KSHamOrb[index] - sePtr[index];
						}
			linalg.invert_square_matrix(tmp, dummy ,false);

			auto gPtr = g.write_phs_grid_ptr_block(ik,iw);
			std::copy(tmp.begin(),tmp.end(),gPtr);
		}
	}
}

} /* namespace gw_flex */
} /* namespace scallop */
