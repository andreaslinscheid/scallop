/*	This file ObservableStatistics.hpp is part of scallop.
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
 *  Created on: Dec 20, 2016
 *      Author: A. Linscheid
 */

#include "scallop/output/ObservableStatistics.h"
#include "scallop/auxillary/BasicFunctions.h"
#include <string>

namespace scallop
{
namespace output
{

template<typename T>
void ObservableStatistics<T>::print_statistics( output::TerminalOut & msg,
		gw_flex::SelfEnergy<T> const& Sigma,
		gw_flex::KohnShamBandStructure<T> const& KS,
		gw_flex::GreensFunctionOrbital<T> const& G,
		bT beta) const
{
	size_t nM = Sigma.get_num_time();

	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	msg << "\n\tSummary of observables:";
	assert( not G.is_in_k_space() );
	assert( not G.is_in_time_space() );
	assert( Sigma.is_in_k_space() );
	assert( not Sigma.is_in_time_space() );

//	msg << "\t(THIS HAS TO BE CORRECTED)Filling per orbital:(l=0..."+std::to_string( G.get_nOrb() - 1)+")";

	bT totalN = bT(0);
	std::string line;
	if ( G.get_spaceGrid_proc().get_proc_index( false, /*iR=*/0 ) == mpi.get_mpi_me())
		for ( size_t l = 0 ; l < G.get_nOrb(); ++l)
		{
			T sum = T(0);
			for (size_t iw =0; iw < G.get_num_time(); ++iw)
				for ( size_t a1 = 0 ; a1 < 2; ++a1)
					for ( size_t s1 = 0 ; s1 < 2; ++s1)
						for ( size_t a2 = 0 ; a2 < 2; ++a2)
							for ( size_t s2 = 0 ; s2 < 2; ++s2)
							{
								sum += auxillary::BasicFunctions::get_v<T>(0)(a1,s1,a2,s2)*G(/*iR=*/0,iw,l,a1,s1,l,a2,s2);
							}
			sum *= 0.5/beta;
			totalN += std::real(sum);
			line += std::string("\t")+std::to_string( std::real(sum) );
		}
	msg << line;
//	msg << "\tTotal filling " << totalN;

//	msg << "\tMagnetiszation per orbital: (l=0..."+std::to_string( G.get_nOrb() - 1)+")";
	auto mtranslate = [] (size_t i) { return i == 1 ? "m_x" : i == 2 ? "m_y" : i == 3 ? "m_z" : "";	};

	bT totalM = bT(0);
	for ( size_t i = 1; i < 4; ++i)
	{
		line = std::string("\t\t")+mtranslate(i)+" per orbital:\n\t";
		if ( G.get_spaceGrid_proc().get_proc_index( false, /*iR=*/0 ) == mpi.get_mpi_me())
			for ( size_t l = 0 ; l < G.get_nOrb(); ++l)
			{
				T sum = T(0);
				for (size_t iw =0; iw < G.get_num_time(); ++iw)
					for ( size_t a1 = 0 ; a1 < 2; ++a1)
						for ( size_t s1 = 0 ; s1 < 2; ++s1)
							for ( size_t a2 = 0 ; a2 < 2; ++a2)
								for ( size_t s2 = 0 ; s2 < 2; ++s2)
								{
									sum += auxillary::BasicFunctions::get_v<T>(i)(a1,s1,a2,s2)*G(/*iR=*/0,iw,l,a1,s1,l,a2,s2);
								}
				sum *= 0.5/beta;
				totalM += std::real(sum);
				line += std::string("\t")+std::to_string( std::real(sum) );
			}
//		msg << line;
	}
//	msg << "\tTotal magnetization " << totalM;

	msg << "\tEstimating mass enhancement factor per orbital (lowest Matsubara frequency of Z):";
	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector V;
	typedef typename auxillary::TemplateTypedefs<bT>::scallop_vector VbT;
	size_t nO = KS.get_nOrb();
	size_t nK = Sigma.get_spaceGrid_proc().get_num_k_grid();
	auxillary::LinearAlgebraInterface<T> linalg;
	gw_flex::MemoryLayout mem;
	mem.initialize_layout_2pt_obj(nO);
	V buffer(nO*2*nO*2);
	VbT ev(nO*2);

	// this stores for each k pt and orbital, the projections onto the
	// input orbitals and the frequency renormalization Z(k,l).
	std::vector<VbT> kResolvedZ0( nO, VbT(nK, bT(0)) );
	for ( size_t ik = 0; ik < nK ; ++ik)
	{
		auto lowestMFreqs = Sigma.read_data_ptr_block(ik,0);
		for ( size_t l1 = 0 ; l1 < nO; ++l1)
			for ( size_t a1 = 0 ; a1 < 2; ++a1)
				for ( size_t s1 = 0 ; s1 < 2; ++s1)
					kResolvedZ0[l1][ik] += 0.25*std::imag(lowestMFreqs[mem.memory_layout_2pt_obj(l1,a1,s1,l1,a1,s1)]);
	}

	for ( size_t io = 0; io < nO ; ++io)
	{
		auto itkm = std::max_element(kResolvedZ0[io].begin(),kResolvedZ0[io].end());
		bT Z = bT(1.0)+ (*itkm)	/ std::imag(auxillary::BasicFunctions::matzubara_frequency_of_index(0,nM,beta));

		msg << "\tMax Z("<<io<<",ik) = "<<Z;
	}
}

} /* namespace output */
} /* namespace scallop */
