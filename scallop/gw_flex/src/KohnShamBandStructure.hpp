/*	This file KohnShamBandStructure.hpp is part of scallop.
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
 *  Created on: Nov 25, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/KohnShamBandStructure.h"
#include "scallop/gw_flex/WannierHamiltonian.h"
#include "scallop/parallel/GridDistribution.h"
#include "scallop/error_handling/Warning.h"
#include "scallop/auxillary/BasicFunctions.h"
#include <memory>
#include <algorithm>

namespace scallop
{
namespace gw_flex
{

template<typename T>
KohnShamBandStructure<T>::KohnShamBandStructure()
{
}

template<typename T>
void KohnShamBandStructure<T>::initialize_from_file(
		std::vector<size_t> grid,
		std::string const & fileName )
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	gridDistr_.distribute_grid(grid);
	size_t nK = gridDistr_.get_num_k_grid();
	typename auxillary::TemplateTypedefs<bT>::scallop_vector kpoints( grid.size()*nK );
	for ( size_t ik = 0; ik < nK; ++ik )
	{
		auto tuple = gridDistr_.k_conseq_local_to_xyz_total( ik );
		for ( size_t id = 0; id < tuple.size() ; ++id)
			kpoints[ik*tuple.size()+id] = bT(tuple[id])/bT(grid[id]);
	}

	useModel_ = false;
	if ( mpi.ioproc() )
	{
		//Read the header and see how we acquire the data for the band structure
		std::ifstream file;
		file.open( fileName.c_str() );
		if ( ! file.good() )
			error_handling::Error(std::string()+"Unable to open file "+fileName,1);
		std::string headerLine;
		file >> headerLine;
		if ( headerLine.compare("model") == 0 )
		{
			useModel_ = true;
		}
		file.close();
	}

	mpi.bcast(useModel_, mpi.ioproc_index() );

	if ( not useModel_ )
	{
		wanHam_.load_wan_ham( fileName );
	}
	else
	{
		//here, every proc opens the file!
		std::ifstream file;
		file.open( fileName.c_str() );
		if ( ! file.good() )
			error_handling::Error(std::string()+"Unable to open file "+fileName,1);
		this->set_model( file );
	}

	this->initialize_layout_2pt_obj( this->get_nOrb() );

	typename auxillary::TemplateTypedefs<T>::scallop_vector unitary;
	this->compute_at_k(kpoints,nK,enk_,unitary);

	akil_.initialize( std::move(grid), this->get_nOrb(), std::move(unitary) );
}

template<typename T>
parallel::GridDistribution<T> const&
KohnShamBandStructure<T>::get_spaceGrid_proc() const
{
	return gridDistr_;
}

template<typename T>
void KohnShamBandStructure<T>::set_model(std::istream & stream)
{
	std::string buffer;
	stream >> buffer; // model

	stream >> buffer; //which model?
	if ( buffer.compare("TwoBandCosine") == 0 )
		model_ = std::make_shared< TwoBandCosine<T> >();

	if ( buffer.compare("OneBandCosine") == 0 )
		model_ = std::make_shared< OneBandCosine<T> >();

	model_->load_model_parameters( stream );
}

template<typename T>
template<class VbT, class V>
void KohnShamBandStructure<T>::compute_at_k(
		VbT const& kpoints,size_t nK,VbT & enk,V & unitary) const
{
	if ( useModel_ )
	{
		vbt kpoints_cpy(kpoints.begin(),kpoints.end());
		vbt enk_cpy;
		v unitary_cpy;
		model_->compute_at_k(kpoints_cpy,nK,unitary_cpy,enk_cpy);
		enk = VbT(enk_cpy.begin(),enk_cpy.end());
		unitary = V(unitary_cpy.begin(),unitary_cpy.end());
	}
	else
	{
		wanHam_.compute_at_k(kpoints,nK,unitary,enk);
	}

	//apply the internal chemical potential offset
	for ( size_t ik = 0 ; ik < nK ; ++ik )
		for ( size_t ib = 0 ; ib < this->get_nOrb() ; ++ib )
			for ( size_t a = 0 ; a < 2 ; ++a )
				for ( size_t is = 0 ; is < 2 ; ++is )
				{
					size_t dIndex = ik*this->get_nOrb()*4 + this->memory_layout_2pt_diagonal(ib,a,is);
					enk[dIndex] -= bT(a == 0 ? 1.0 :-1.0)*mu_;
				}

	assert( enk.size() == nK*this->get_nOrb()*4 );
	assert( unitary.size() == nK*this->get_nOrb()*4*this->get_nOrb()*4 );
}

template<typename T>
size_t KohnShamBandStructure<T>::get_nOrb() const
{
	if ( useModel_ )
		return model_->get_nOrb();
	return wanHam_.get_nOrb();
}

template<typename T>
typename KohnShamBandStructure<T>::bT
KohnShamBandStructure<T>::operator ()(size_t ik, size_t n, size_t spinCnl, size_t phCnl) const
{
	size_t dIndex = ik*this->get_nOrb()*4 + this->memory_layout_2pt_diagonal(n,phCnl,spinCnl);
	ASSERT( dIndex < enk_.size(),\
			std::string("Out of range: index is ")\
			+std::to_string(dIndex)+" while range is "+std::to_string(enk_.size()) );
	return enk_[dIndex];
}

template<typename T>
UnitaryWannierKSBands<T> const &
KohnShamBandStructure<T>::get_unitary() const
{
	return akil_;
}

template<typename T>
typename KohnShamBandStructure<T>::vbt const &
KohnShamBandStructure<T>::get_bands() const
{
	return enk_;
}

template<typename T>
void KohnShamBandStructure<T>::set_chem_pot( bT chemicalPotential )
{
	bT shift = (chemicalPotential-mu_);
	size_t nK = gridDistr_.get_num_k_grid();
	size_t nB = this->get_nOrb();

	for ( size_t ik = 0 ; ik < nK ; ++ik )
		for ( size_t ib = 0 ; ib < nB ; ++ib )
			for ( size_t a = 0 ; a < 2 ; ++a )
				for ( size_t is = 0 ; is < 2 ; ++is )
				{
					size_t dIndex = ik*this->get_nOrb()*4 + this->memory_layout_2pt_diagonal(ib,a,is);
					enk_[dIndex] -= bT(a == 0 ? 1.0 :-1.0)*shift;
				}
	mu_ = chemicalPotential;
}

template<typename T>
typename KohnShamBandStructure<T>::bT
KohnShamBandStructure<T>::get_chem_pot( ) const
{
	return mu_;
}

template<typename T>
void KohnShamBandStructure<T>::adjust_filling( bT numElectrons, bT invTemp )
{
	auxillary::BasicFunctions bf;
	auto oldMu = mu_;

	struct Filling
	{
		Filling( KohnShamBandStructure<T> * ptr, bT invTemp, bT Ne) :
			Ne_(Ne), ptr_(ptr) , invTemp_(invTemp){};

		bT operator() (bT x) const
		{
			ptr_->set_chem_pot( x );
			bT nElec = ptr_->compute_N_electrons( invTemp_ );
			return (nElec - Ne_);
		};

		bool check_conv(bT f_xn, bT f_xnp1) const
		{
			return std::abs(f_xn - f_xnp1) < 0.0000001;
		};

	private:
		bT Ne_;
		KohnShamBandStructure<T> * ptr_;
		bT invTemp_ ;
	};

	Filling filling( this, invTemp, numElectrons );

	auto energyRange = this->get_band_width();

	std::pair<bool,bT> result;
	bf.secant_method( bT((energyRange.first+energyRange.second)/2.0),
			energyRange,
			filling,1000,
			result);

	if ( result.first )
	{
		this->set_chem_pot( result.second );
	}
	else
	{
		error_handling::Warning warn( std::string("Cannot find chemical potential to yield ")
			+std::to_string( numElectrons ) + "electrons. Using previous chemical potential ...");
		this->set_chem_pot( oldMu );
	}
}

template<typename T>
std::pair<typename KohnShamBandStructure<T>::bT,typename KohnShamBandStructure<T>::bT>
KohnShamBandStructure<T>::get_band_width() const
{
	auto result = std::make_pair( bT(0), bT(0) );
	decltype(enk_) nambu11Part( enk_.size()/2 );
	size_t nK = gridDistr_.get_num_k_grid();
	size_t nB = this->get_nOrb();
	size_t c =0;
	for ( size_t ik = 0 ; ik < nK ; ++ik )
		for ( size_t ib = 0 ; ib < nB ; ++ib )
			for ( size_t is = 0 ; is < 2 ; ++is )
			{
				//We use only the Nambu 11 compontent
				size_t dIndex = ik*this->get_nOrb()*4 + this->memory_layout_2pt_diagonal(ib,0,is);
				nambu11Part[c++] = enk_[dIndex];
			}
	auto pairit = std::minmax_element( nambu11Part.begin(), nambu11Part.end() );
	if( (!(pairit.first == nambu11Part.end())) && (!(pairit.second == nambu11Part.end())) )
		result = std::make_pair(*pairit.first,*pairit.second);
	return result;
}

template<typename T>
typename KohnShamBandStructure<T>::bT
KohnShamBandStructure<T>::compute_N_electrons(bT invTemperature, bT add_shift_chem_pot) const
{
	bT result = 0;
	size_t nK = gridDistr_.get_num_k_grid();
	size_t nB = this->get_nOrb();
	for ( size_t ik = 0 ; ik < nK ; ++ik )
		for ( size_t ib = 0 ; ib < nB ; ++ib )
			for ( size_t is = 0 ; is < 2 ; ++is )
			{
				//We use only the Nambu 11 compontent
				size_t dIndex = ik*this->get_nOrb()*4 + this->memory_layout_2pt_diagonal(ib,0,is);
				bT e = enk_[dIndex] - add_shift_chem_pot;
				result += auxillary::BasicFunctions::fermi_funcs(invTemperature,e);
			}
	result /= gridDistr_.get_num_grid();

	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	mpi.sum(result);
	return result;
}

template<typename T>
void KohnShamBandStructure<T>::get_local_orbital_KSHam(size_t ik,
		v & orbitalHamiltonian,
		auxillary::LinearAlgebraInterface<T> const& linalg
		)const
{
	v buffer1;
	v buffer2;
	this->get_local_orbital_KSHam(ik,orbitalHamiltonian,linalg,buffer1,buffer2);
}

template<typename T>
void KohnShamBandStructure<T>::get_local_orbital_KSHam(size_t ik,
		v & orbitalHamiltonian,
		v & buffer1,
		v & buffer2,
		auxillary::LinearAlgebraInterface<T> const& linalg) const
{
	size_t nB = std::pow(this->get_nOrb()*4,2);

	if ( orbitalHamiltonian.size() != nB)
		orbitalHamiltonian = v(nB);

	if ( buffer1.size() != nB)
		buffer1 = v(nB);

	if ( buffer2.size() != nB)
		buffer2 = v(nB);

	//Construct the KS bands in Orbital space at this k pt
	for ( size_t m1 = 0 ; m1 < nB ; ++m1)
		for ( size_t m2 = 0 ; m2 < nB ; ++m2)
			buffer1[m1*nB+m2] = std::conj(akil_(ik,m1,m2));

	linalg.matrix_times_diagonal_matrix( akil_.read_phs_grid_ptr_block(ik), nB, &( enk_[ik*nB] ), buffer2.data() );
	linalg.matrix_times_matrix( buffer1.data(), nB, buffer2.data(), orbitalHamiltonian.data() );
}

} /* namespace gw_flex */
} /* namespace scallop */
