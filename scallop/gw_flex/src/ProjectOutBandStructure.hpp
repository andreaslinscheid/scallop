/*	This file ProjectOutBandStructure.hpp is part of scallop.
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
 *  Created on: Dec 8, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/ProjectOutBandStructure.h"
#include <algorithm>

namespace scallop
{
namespace gw_flex
{
template<typename T>
ProjectOutBandStructure<T>::ProjectOutBandStructure(
		KohnShamBandStructure<T> bands) : bands_(std::move(bands))
{

}

template<typename T>
size_t ProjectOutBandStructure<T>::size_per_block() const
{
	return 16*bands_.get_nOrb()*bands_.get_nOrb();
}

template<typename T>
size_t ProjectOutBandStructure<T>::size_per_block_after() const
{
	return 1;
}

template<typename T>
void ProjectOutBandStructure<T>::project_out_KS_bands_half_freq(
		SelfEnergy<T> const& SE,
		V const& MatsFreq,
		V & mappedMatsFreq, container_type & mappedData) const
{
	assert( ! SE.is_in_time_space() );
	assert( SE.is_in_k_space() );
	assert( MatsFreq.size() == SE.get_num_time() );
	assert( SE.get_data_block_size() == this->size_per_block() );

	size_t nK = SE.get_spaceGrid_proc().get_num_k_grid();
	size_t nM = SE.get_num_time();
	size_t nB = this->size_per_block();

	if ( mappedMatsFreq.size() != nM/2 )
		mappedMatsFreq = V(nM/2);
	std::copy(MatsFreq.begin()+nM/2,MatsFreq.end(),mappedMatsFreq.begin());

	if ( mappedData.size() != nK*nM/2*nB )
		mappedData = container_type(nK*nM/2*nB);

	auto unitary = bands_.get_unitary();
	assert( unitary.get_spaceGrid_proc().get_grid() == SE.get_spaceGrid_proc().get_grid());

	//for each point, copy the second half of the data array corresponding to the
	// negative frequencies
	for ( size_t ik = 0 ; ik <nK; ++ik )
	{
		T const * se_ptr = SE.read_phs_grid_ptr_block(ik,nM/2);
		std::copy(se_ptr,se_ptr+nM/2*nB, mappedData.begin()+ik*nM/2*nB );
	}
}

template<typename T>
void ProjectOutBandStructure<T>::before(container_type & selfEnergyAlongPath,
		parallel::IrregularGridDistribution<T> const& path,
		V const& MatsGrid) const
{

}

template<typename T>
void ProjectOutBandStructure<T>::after( container_type & dataPath,
		parallel::IrregularGridDistribution<T> const& path,
		V const& omega) const
{
	//To extract the particle band structure we have to compute the antihermitian
	//	part of the Nambu 1,1 component and diagonalize it.
	size_t nB = this->size_per_block();
	size_t nO = bands_.get_nOrb();
	size_t nOs = bands_.get_nOrb()*2;
	size_t nOns = bands_.get_nOrb()*4;
	size_t nK = path.get_n_pts();
	assert( nOs*nOs*4 == nB);

	//local typdefs derived from the contain type
	typedef typename container_type::value_type cT;
	typedef typename auxillary::TypeMapComplex<cT>::type cbT;
	typedef typename auxillary::TemplateTypedefs<cbT>::scallop_vector cVbT;

	MemoryLayout meml;
	meml.initialize_layout_2pt_obj(nO);

	//obtain the k points, then the KS eigenvalues and the unitary matrices
	size_t dim = path.get_vector( 0 ).size();
	cVbT kpoints( nK * dim);
	for (size_t ik = 0 ; ik < nK; ++ik)
	{
		auto k = path.get_vector( ik );
		for (size_t i = 0 ; i < dim ;++i)
			kpoints[ik*dim+i] = k[i];
	}
	cVbT enk;
	container_type unitaryPath;
	bands_.compute_at_k(kpoints,nK,enk,unitaryPath);

	auxillary::LinearAlgebraInterface<cT> linalg;

	container_type localCpyNambu11(nOs*nOs);
	container_type KSOrb(nB);
	container_type conjUnitary(nB);
	container_type bnd(nOns);
	container_type tmp(nB);
	cbT dummy;

	container_type localCpy(nB);
	cVbT ev(nOs);

	container_type mappingData( nK*omega.size() );

	for ( size_t ik = 0 ; ik <nK; ++ik )
	{
		//Construct the KS bands in Orbital space at this k pt
		for ( size_t m1 = 0 ; m1 < nOns ; ++m1)
			for ( size_t m2 = 0 ; m2 < nOns ; ++m2)
				conjUnitary[m1*nOns+m2] = std::conj(unitaryPath[(ik*nOns+m2)*nOns+m1]);

		for ( size_t n = 0 ; n < nOns; ++n)
			bnd[n] = enk[ik*nOns+n];
		auto unitaryThisK = &( unitaryPath[ik*nOns*nOns] );
		linalg.matrix_times_diagonal_matrix( unitaryThisK, nOns, bnd.data(), tmp.data() );
		linalg.matrix_times_matrix( conjUnitary.data(), nOns, tmp.data(), KSOrb.data() );

		for ( size_t iw = 0 ; iw < omega.size(); ++iw )
		{
			//Construct the inverse GF on the nearly real axis
			auto it = dataPath.begin()+(ik*omega.size()+iw)*nB;
			std::copy( it, it+nB, localCpy.begin() );

			for ( size_t i = 0 ; i < nOns; ++i)
				for ( size_t j = 0 ;j < nOns; ++j)
					localCpy[i*nOns+j] = cT((i==j?omega[iw]:T(0)))-KSOrb[i*nOns+j]-localCpy[i*nOns+j];

			//invert to obtain the GF
			linalg.invert_square_matrix(localCpy,dummy);

			std::copy( localCpy.begin(),localCpy.end(), it );

			//Obtain the antihermititan part and eigenvalues
			size_t a1=0,a2=0;
			for ( size_t l1 = 0 ; l1 < nO; ++l1 )
				for ( size_t s1 = 0 ; s1 < 2; ++s1 )
					for ( size_t l2 = 0 ; l2 < nO; ++l2 )
						for ( size_t s2 = 0 ; s2 < 2; ++s2 )
							localCpyNambu11[((l1*2+s1)*nO+l2)*2+s2] = *(it+meml.memory_layout_2pt_obj(l1,a1,s1,l2,a2,s2));

			for ( size_t i = 0 ; i < nOs; ++i)
				for ( size_t j = i ; j < nOs; ++j)
				{
					cT antiHermititan = cbT(0.5)*(localCpyNambu11[i*nOs+j] - std::conj(localCpyNambu11[j*nOs+i]));
					localCpyNambu11[i*nOs+j] = antiHermititan;
					localCpyNambu11[j*nOs+i] = antiHermititan;
				}

			//If A is skew Hermitian, then -i*A is hermititan
			for ( auto &e : localCpyNambu11 )
				e *= -1.0*T(0,1);

			linalg.hermitian_eigensystem(true,false,localCpyNambu11.data(),nOs,ev.data());

			cT sum_ev=cT(0);
			for(auto e: ev)
				sum_ev += cT(0,1)*e;
			mappingData[ik*omega.size()+iw] = sum_ev;
		}
	}
	dataPath=std::move(mappingData);
}

} /* namespace gw_flex */
} /* namespace scallop */
