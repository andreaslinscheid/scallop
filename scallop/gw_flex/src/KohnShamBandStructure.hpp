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
#include <memory>

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
	for ( auto &e : enk )
		e -= mu_;

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
	for ( auto &e : enk_ )
		e -= shift;
	mu_ = chemicalPotential;
}

template<typename T>
void KohnShamBandStructure<T>::adjust_filling( bT numElectronPerSpin )
{

}

} /* namespace gw_flex */
} /* namespace scallop */
