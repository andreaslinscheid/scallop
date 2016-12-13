/*	This file WannierHamiltonian.hpp is part of scallop.
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
 *  Created on: Nov 26, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/WannierHamiltonian.h"
#include "scallop/parallel/MPIModule.h"
#include "scallop/output/TerminalOut.h"
#include "scallop/auxillary/LinearAlgebraInterface.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
void WannierHamiltonian<T>::load_wan_ham( std::string const & fileNameWannierHam )
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	output::TerminalOut msg;

	if ( mpi.ioproc() )
	{
		//format is described in the Wannier90 manual
		std::ifstream file;
		file.open( fileNameWannierHam.c_str() );
		if ( ! file.good() )
			error_handling::Error(std::string()+"Unable to open file "+fileNameWannierHam,1);
		msg << "Reading file " << fileNameWannierHam;

		//Header, we discard the first line that says when the file was created
		std::string s;
		std::getline(file,s);
		file >> nOrb_;
		size_t nRVects;
		file >> nRVects;
		msg << "\tThe number of orbitals is " << nOrb_
				<< "; found " << nRVects << " R vectors in the file.";

		weights_ = vbT(nRVects, bT(0) );
		RGrid_ = vbT(3*nRVects, bT(0) );
		wanHam_ = v( nOrb_*nOrb_*nRVects );

		for ( size_t i = 0 ; i < nRVects; ++i )
			file >> weights_[i];

		std::vector<bT> R(3);
		int mu,nu;
		bT hopRe, hopImag;
		auto set_wanHam_from_line = [&] ()
			{
				bool read_success= true;
				for ( auto &&rxi : R )
					read_success = read_success && (file >> rxi);
				read_success = read_success && (file >> mu >> nu >> hopRe >> hopImag) ;
				//nu and mu are in the fortran convenction 1...nOrb
				return read_success;
			};

		for ( size_t iR = 0 ; iR < nRVects; ++iR)
		{
			for ( size_t i = 0 ; i < nOrb_*nOrb_; ++i)
			{
				if ( ! set_wanHam_from_line() )
					error_handling::Error("Wannier Hamiltonian data file "+fileNameWannierHam
							+" could not be read successfully!\n Failed at iR="+std::to_string(iR)
							+", mu="+std::to_string(mu)+", nu="+std::to_string(nu)
							+", i.e. i="+std::to_string(i));
				wanHam_[(iR*nOrb_+(nu-1))*nOrb_+(mu-1)] = T(hopRe,hopImag);
			}
			std::copy(R.begin(),R.end(),RGrid_.begin()+iR*3);
		}
	}

	mpi.bcast(weights_, mpi.ioproc_index() );
	mpi.bcast(RGrid_, mpi.ioproc_index() );
	mpi.bcast(wanHam_, mpi.ioproc_index() );
	mpi.bcast(nOrb_, mpi.ioproc_index() );
};

template<typename T>
template<class VbT, class V>
void WannierHamiltonian<T>::compute_at_k( VbT kpts, size_t nkpts, V & unitary, VbT & energyEV ) const
{
	size_t dim = kpts.size()/nkpts;
	unitary = V( 16*nkpts*nOrb_*nOrb_, typename V::value_type(0) );
	energyEV = VbT( 4*nkpts*nOrb_ );

	auxillary::TemplateTypedefs< std::complex<double> >::scallop_vector hamltonianAtK( 16*nkpts*nOrb_*nOrb_, std::complex<double>(0) );
	for ( size_t ik = 0 ; ik < nkpts; ++ik)
	{
		for ( size_t iR = 0 ; iR < weights_.size() ; ++iR)
		{
			double dprod = 0.0;
			for ( size_t id = 0 ; id < dim; ++id )
				dprod += RGrid_[iR*3+id]*kpts[ik*dim+id];
			std::complex<double> phase = std::exp( std::complex<double>( 0, 2*M_PI*dprod ) );
			for ( size_t mu = 0 ; mu < nOrb_ ; ++mu)
				for ( size_t nu = 0 ; nu < nOrb_ ; ++nu)
					hamltonianAtK[(ik*nOrb_+mu)*nOrb_+nu] += wanHam_[(iR*nOrb_+nu)*nOrb_+mu]/weights_[iR]*phase;
		}
	}

	MemoryLayout meml;
	meml.initialize_layout_2pt_obj( nOrb_ );

	auxillary::LinearAlgebraInterface<std::complex<double> > linAlgebra;
	auxillary::TemplateTypedefs<double>::scallop_vector ev( nOrb_ );
	for ( size_t ik = 0 ; ik < nkpts; ++ik)
	{
		linAlgebra.hermitian_eigensystem(true,true,&(hamltonianAtK[ik*nOrb_*nOrb_]),nOrb_,ev.data());
		//From here on kLocalHamiltonian contains the eigenvectors

		for ( size_t mu = 0 ; mu < nOrb_; ++mu)
			for ( size_t a = 0 ; a < 2; ++a)
				for ( size_t s = 0 ; s < 2; ++s)
				energyEV[ik*4*nOrb_ + meml.memory_layout_2pt_diagonal(mu,a,s)]
				         = 1000.0*//eV to meV
				         	 (a == 0 ? 1.0 : -1.0 )*static_cast<bT>(ev[mu]);

		auto it = &( unitary[ik*nOrb_*nOrb_*16] );
		std::fill(it,it+nOrb_*nOrb_*16,typename V::value_type(0));
		//Lapack return the eigenvectors in Column major.
		for ( size_t mu = 0 ; mu < nOrb_; ++mu)
			for ( size_t ns1 = 0 ; ns1 < 4; ++ns1)
				for ( size_t nu = 0 ; nu < nOrb_; ++nu)
					unitary[ik*nOrb_*nOrb_*16+meml.memory_layout_2pt_obj_nsc(mu,ns1,nu,ns1)] = hamltonianAtK[(ik*nOrb_+nu)*nOrb_+mu];
	}
}

template<typename T>
size_t WannierHamiltonian<T>::get_nOrb() const
{
	return nOrb_;
}


} /* namespace gw_flex */
} /* namespace scallop */
