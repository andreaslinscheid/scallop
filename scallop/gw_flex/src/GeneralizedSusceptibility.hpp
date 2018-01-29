/*	This file GeneralizedSusceptibility.hpp is part of scallop.
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
 *  Created on: Nov 11, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/GeneralizedSusceptibility.h"
#include <assert.h>

namespace scallop
{
namespace gw_flex
{

template<typename T>
GeneralizedSusceptibility<T>::GeneralizedSusceptibility() :
	MatsubaraImagTimeFourierTransform<T>( /* bool Fermi=*/ false )
{

}

template<typename T>
void GeneralizedSusceptibility<T>::initialize_zero(
		size_t nM,
		size_t nO,
		size_t channels,
		std::vector<size_t> grid,
		bool intimeSpace,
		bool inKSpace)
{
	MemoryLayout::initialize_layout_4pt_scalar_obj( nO ,channels);
	size_t nC = this->get_nChnls();

	size_t blocksize = std::pow( nO , 4 )*nC*nC;
	typename auxillary::TemplateTypedefs<T>::scallop_vector data;

	// initialize an object with data being zero
	this->initialize( nM, grid, blocksize, intimeSpace, inKSpace, data );
	bufferQ1_ = bufferQ2_ = typename auxillary::TemplateTypedefs<T>::scallop_vector( nC*nO*nO* 16 );
	bufferSBlock_ = typename auxillary::TemplateTypedefs<T>::scallop_vector( nC*nC*std::pow(nO,4) );
}

template<typename T>
void GeneralizedSusceptibility<T>::plot_static(
		output::SusceptibilityPlotter & plotter) const
{
	assert( this->is_in_k_space() );
	assert( not this->is_in_time_space() );

	size_t channels =  this->get_nChnls();
	size_t nK = this->get_spaceGrid_proc().get_num_k_grid();
	int blockDim = static_cast<int>( std::pow(this->get_nOrb(),2) * channels );

	bool plotSpin = (channels > 1) or (channelsOffset_ > 0);
	bool plotCharge = not plotSpin;
	if ( (plotSpin and plotter.do_plot_spin() ) or (plotCharge and plotter.do_plot_charge() ) )
	{
		plotter.set_buffer( this->get_spaceGrid_proc().get_k_grid(),
				this->get_nOrb(),blockDim*blockDim*nK,channels);
	}

	for ( size_t ik = 0 ; ik < nK; ++ik)
	{
		auto susct_d_begin = this->read_phs_grid_ptr_block(ik,/*iw=*/0);
		auto susct_d_end = this->read_phs_grid_ptr_block(ik,/*iw=*/1);
		typename auxillary::TemplateTypedefs<T>::scallop_vector data(susct_d_begin,susct_d_end);
		plotter.append_data_to_buffer(data);
	}

	if ( (plotSpin and plotter.do_plot_spin() ) or (plotCharge and plotter.do_plot_charge() ) )
	{
		plotter.plot_static_k();
	}
}

template<typename T>
void GeneralizedSusceptibility<T>::compute_from_gf(
		GreensFunctionOrbital<T> const & GF,
		size_t channels,
		size_t channel_offset)
{
	//check if we need to update, or initialize
	if ( ! this->is_init() )
		this->initialize_zero(
				GF.get_num_time(), GF.get_nOrb(), channels, GF.get_spaceGrid_proc().get_grid(),
				/*In time space=*/ true, /*In k space=*/ false );

	channelsOffset_ = channel_offset;

	size_t nO = this->get_nOrb();
	size_t nC = this->get_nChnls();

	assert( (not this->is_in_k_space()) and this->is_in_time_space() );
	assert( nO ==  GF.get_nOrb() );
	assert( this->get_num_time() ==  GF.get_num_time() );
	assert( this->get_spaceGrid_proc().get_grid() ==  GF.get_spaceGrid_proc().get_grid() );
	assert( (not GF.is_in_k_space()) and GF.is_in_time_space() );

	MemoryLayout gf_layout;
	gf_layout.initialize_layout_2pt_obj( nO );

	/*
	 * NOTE: (TODO)
	 *
	 * We do not implement non-inversionally symmetric systems at this point.
	 * This would open a can of worms, since the inverse index may well be on another processor.
	 * Similarly, for non-local interactions, where this routine would have to be called with an R-offset.
	 * If both is not the case, then this is a GF local product (in space)!
	 */
	for ( size_t iR = 0; iR < this->get_spaceGrid_proc().get_num_R_grid(); ++iR)
	{
		for ( size_t it = 0; it < this->get_num_time(); ++it)
		{
			//Note that the shift here (-\tau) -> (beta - \tau) causes a minus sign to appear because
			// of the Fermionic Matsubara rules.
			size_t it_inverse = static_cast<size_t>(static_cast<int>( this->get_num_time() ) - static_cast<int>(it) - 1 );

			T * blockPtrSusct = this->write_phs_grid_ptr_block(iR,it);
			T const * blockPtrG1 = GF.read_phs_grid_ptr_block(iR,it_inverse);
			T const * blockPtrG2 = GF.read_phs_grid_ptr_block(iR,it);

			//Construct the 4Norbital x 16 matrix Q1 and the 16 x 4Norbital matrix Q2
			//We use a hard-coded matrix multiplication v_j * G
			for ( size_t j = 0 ; j < nC; ++j )
				for ( size_t l1 = 0 ; l1 < nO; ++l1 )
					for ( size_t a1 = 0 ; a1 < 2; ++a1 )
						for ( size_t s1 = 0 ; s1 < 2; ++s1 )
							for ( size_t l2 = 0 ; l2 < nO; ++l2 )
								for ( size_t a2 = 0 ; a2 < 2; ++a2 )
									for ( size_t s2 = 0 ; s2 < 2; ++s2 )
									{
										T prefactorG1;
										size_t gf1_index_NS_transpose = 0;
										this->v_matrix_multiplication( j+channelsOffset_, l1, a2, s2, l2, a1, s1,
												gf_layout,gf1_index_NS_transpose ,prefactorG1);

										size_t Q1_index = (((((j*nO+l1)*nO+l2)*2+a1)*2+s1)*2+a2)*2+s2;
										assert( Q1_index < bufferQ1_.size() );
										bufferQ1_[Q1_index] = prefactorG1*blockPtrG1[ gf1_index_NS_transpose ];

										T prefactorG2;
										size_t gf2_index_NS = 0;
										this->v_matrix_multiplication( j+channelsOffset_, l1, a1, s1, l2, a2, s2,
												gf_layout,gf2_index_NS ,prefactorG2);

										size_t Q2_index = (((((a1*2+s1)*2+a2)*2+s2)*nC+j)*nO+l1)*nO+l2;
										assert( Q2_index < bufferQ2_.size() );
										bufferQ2_[Q2_index] =  prefactorG2*blockPtrG2[ gf2_index_NS ];
									}

			int dimSust = static_cast<int>(nC*nO*nO);
			linAlgModule_.call_gemm(false, false,
					dimSust,dimSust,16,
					T(-1.0/4.0),	bufferQ1_.data(), 16,
								bufferQ2_.data(), dimSust,
		            T(0.0),bufferSBlock_.data(),dimSust);

			//Now, the susceptibility contains the right data, but not in the right order.
			//The current order is j,l1,l3,jp,l4,l2 if the first GF has orbital indices l1,l3 and the second l4,l2
			//We rearrange this to j,jp,l1,l2,l3,l4 in the previous nomenclature
			for ( size_t j = 0 ; j < nC; ++j )
				for ( size_t l1 = 0 ; l1 < nO; ++l1 )
					for ( size_t l3 = 0 ; l3 < nO; ++l3 )
						for ( size_t jp = 0 ; jp < nC; ++jp )
							for ( size_t l4 = 0 ; l4 < nO; ++l4 )
								for ( size_t l2 = 0 ; l2 < nO; ++l2 )
								{
									size_t bufferIndex = ((((j*nO+l1)*nO+l3)*nC+jp)*nO+l4)*nO+l2;
									size_t susctIndex = ((((j*nC+jp)*nO+l1)*nO+l2)*nO+l3)*nO+l4;
									blockPtrSusct[susctIndex] = bufferSBlock_[bufferIndex];
								}
		}
	}
}

template<typename T>
void GeneralizedSusceptibility<T>::v_matrix_multiplication(
		size_t j, size_t l1, size_t a1, size_t s1, size_t l2,  size_t a2, size_t s2,
		MemoryLayout const& gf_layout,
		size_t &gf_index, T & prefactor) const
{
	auto flip_spin = [] (size_t s) { return s == 0 ? 1 : 0; };

	switch (j) {
		case 0:
			prefactor = a1 != 0 ? T(1) : T(-1);
			gf_index = gf_layout.memory_layout_2pt_obj(l1,a1,s1,l2,a2,s2);
			break;
		case 1:
			prefactor = a1 != 0 ? T(1) : T(-1);
			gf_index = gf_layout.memory_layout_2pt_obj(l1,a1,flip_spin(s1),l2,a2,s2);
			break;
		case 2:
			prefactor = s1 != 0 ? T(0,1) : T(0,-1);
			gf_index = gf_layout.memory_layout_2pt_obj(l1,a1,flip_spin(s1),l2,a2,s2);
			break;
		case 3:
			prefactor = (s1 + a1 == 1) ? T(-1) : T(1);
			gf_index = gf_layout.memory_layout_2pt_obj(l1,a1,s1,l2,a2,s2);
			break;
		default:
			break;
	}
}

template<typename T>
void GeneralizedSusceptibility<T>::RPA_enhancement(
		InteractionMatrix<T> const& interMat)
{
	AdiabaticUpscale a;
	this->RPA_enhancement(interMat,a);
}

template<typename T>
void GeneralizedSusceptibility<T>::RPA_enhancement(
		InteractionMatrix<T> const& interMat,
		AdiabaticUpscale & a )
{
	output::SusceptibilityPlotter p;
	this->RPA_enhancement(interMat,a,p);
}

template<typename T>
void GeneralizedSusceptibility<T>::RPA_enhancement(
		InteractionMatrix<T> const& interMat,
		output::SusceptibilityPlotter & plotter)
{
	AdiabaticUpscale a;
	this->RPA_enhancement(interMat,a,plotter);
}

template<typename T>
void GeneralizedSusceptibility<T>::RPA_enhancement(
		InteractionMatrix<T> const& interMat,
		AdiabaticUpscale & a,
		output::SusceptibilityPlotter & plotter)
{
	assert( this->get_nOrb() == interMat.get_nOrb() );
	assert( this->get_nChnls() == interMat.get_nChnls() );
	assert( this->is_in_k_space() );
	assert( not this->is_in_time_space() );

	size_t channels =  this->get_nChnls();
	size_t nK = this->get_spaceGrid_proc().get_num_k_grid();
	size_t nO = this->get_nOrb();
	int blockDim = static_cast<int>( std::pow(nO,2) * channels );

	bool plotSpin = (channels > 1) or (channelsOffset_ > 0);
	bool plotCharge = not plotSpin;
	bool recordData = false;
	if ( (plotSpin and plotter.do_plot_spin() ) or (plotCharge and plotter.do_plot_charge() ) )
	{
		plotter.set_buffer( this->get_spaceGrid_proc().get_k_grid(),
				this->get_nOrb(),blockDim*blockDim*nK,channels);
		recordData = true;
	}

	typename auxillary::TemplateTypedefs<T>::scallop_vector chiTimesI(blockDim*blockDim,T(0));
	auto IdPlusChiTimesI = chiTimesI;
	auto IdPlusChiTimesIInv = chiTimesI;
	auto interactionScaled = chiTimesI;
	auto enhChiTimesI = chiTimesI;

	typename auxillary::TemplateTypedefs<bT>::scallop_vector eigenvalues;
	bool positiveDefinite;
	bool beyondInstability = false;
	std::vector<size_t> unstableK;
	std::vector<bT> unstableMinEV;

	gw_flex::MemoryLayout ml_inter;
	ml_inter.initialize_layout_4pt_scalar_obj(nO,channels);

	for ( size_t ik = 0 ; ik < nK; ++ik)
	{
		for ( size_t iw = 0 ; iw < this->get_num_time(); ++iw)
		{
			auto susct = this->write_phs_grid_ptr_block(ik,iw);
			for ( size_t j = 0 ; j < channels; ++j )
				for ( size_t jp = 0 ; jp < channels; ++jp )
					for ( size_t l1 = 0 ; l1 < nO; ++l1 )
						for ( size_t l2 = 0 ; l2 < nO; ++l2 )
							for ( size_t l3 = 0 ; l3 < nO; ++l3 )
								for ( size_t l4 = 0 ; l4 < nO; ++l4 )
									interactionScaled[ml_inter.memory_layout_4pt_scalar_obj(j,jp,l1,l2,l3,l4)]
									                  = interMat(j,jp,l1,l2,l3,l4)*a.get_scaling();

			linAlgModule_.call_gemm(false, false,
					blockDim,blockDim,blockDim,
					T(4.0),susct, blockDim,
							interactionScaled.data(), blockDim,
		            T(0.0),chiTimesI.data(),blockDim);

			std::copy(chiTimesI.data(), chiTimesI.data()+chiTimesI.size(), IdPlusChiTimesI.data() );
			for ( int i = 0 ; i < blockDim; ++i)
				for ( int j = 0 ; j < blockDim; ++j)
					IdPlusChiTimesI[i*blockDim+j] = T(i==j?1.0:0.0)-IdPlusChiTimesI[i*blockDim+j];

			this->get_linAlg_module().check_definite_invert_hermitian_matrix(
					IdPlusChiTimesI,IdPlusChiTimesIInv, positiveDefinite, eigenvalues );

			if ( not positiveDefinite )
			{
				beyondInstability = true;
				unstableK.push_back(ik);
				auto itmin = std::min_element(eigenvalues.begin(),eigenvalues.end());
				unstableMinEV.push_back( (*itmin) );
			}

			//Record the data in case we want to plot it
			if ( recordData and (iw == 0) )
			{
				//IdPlusChiTimesI serves as a buffer here that is temporarily overwritten
				//with the RPA susceptibility.
				std::copy(susct,susct+blockDim*blockDim, IdPlusChiTimesI.data() );
				this->get_linAlg_module().matrix_times_matrix(
						IdPlusChiTimesIInv.data(),blockDim,
						enhChiTimesI.data(),
						IdPlusChiTimesI.data() );
				plotter.append_data_to_buffer(IdPlusChiTimesI);
			}

			//double counting correction
			for ( int i = 0 ; i < blockDim; ++i)
				IdPlusChiTimesIInv[i*blockDim+i] -= 0.5;

			this->get_linAlg_module().matrix_times_matrix(
					IdPlusChiTimesIInv.data(),blockDim,
					chiTimesI.data(),
					enhChiTimesI.data() );

			//Currently, we include a factor of two that comes from the ladder series
			linAlgModule_.call_gemm(false, false,
					blockDim,blockDim,blockDim,
					T(-8.0),interactionScaled.data(), blockDim,
							enhChiTimesI.data(), blockDim,
					T(0.0),susct,blockDim);
		}
	}

	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	if ( mpi.false_on_all(  beyondInstability ) )
	{
		output::TerminalOut msg;
		auto itmin = std::min_element(unstableMinEV.begin(),unstableMinEV.end());
		bT minEV = (itmin == unstableMinEV.end() ? 1.0 : (*itmin));
		size_t procMin;
		mpi.min(minEV,procMin);
		std::string type = this->get_nChnls() == 1 ? "charge" : "spin";
		std::string line = std::string("\tInstability in the ")+type+" channel at k pt:\t";
		std::vector<size_t> k;
		if ( procMin == mpi.get_mpi_me() )
		{
			int ikm = itmin - unstableMinEV.begin();
			assert ( (ikm >= 0) && ( static_cast<size_t>(ikm) < unstableK.size() ) );
			k = this->get_spaceGrid_proc().k_conseq_local_to_xyz_total( unstableK[ikm] );
		}
		mpi.bcast(k,procMin);
		auto grid = this->get_spaceGrid_proc().get_grid();
		std::vector<bT> kvector(k.size());
		for ( size_t id = 0 ; id < k.size(); ++id )
		{
			kvector[id] = bT(k[id])/bT(grid[id]);
			line += std::to_string( kvector[id] ) + "\t";
		}
		msg << line << "\n\tSmallest eigenvalue: " << minEV;
		//Currently we don't track the magnetic/orbital eigenvector.
		a.track_instability(kvector,std::vector<T>(eigenvalues.size()));

		if ( a.is_soft() && a.steps_ok() )
		{
			msg << "\tReducing the interaction by a factor " << a.get_scaling();
		}
		else
		{
			msg << "\tNo further reduction of the interaction, we are beyond the instability.";
		}
	}
	else
	{
		if ( a.is_soft() )
			a.set_stable();
	}

	if ( (plotSpin and plotter.do_plot_spin() ) or (plotCharge and plotter.do_plot_charge() ) )
		plotter.plot_static_k();
}

template<typename T>
auxillary::LinearAlgebraInterface<T> const&
GeneralizedSusceptibility<T>::get_linAlg_module() const
{
	return linAlgModule_;
}


template<typename T>
T & GeneralizedSusceptibility<T>::operator() (size_t ik, size_t iw, size_t j, size_t jp, size_t m1,  size_t m2)
{
	return *(this->write_phs_grid_ptr_block(ik,iw)+this->memory_layout_combined_notation_4pt_scalar_obj(j,jp,m1,m2));
}

template<typename T>
T GeneralizedSusceptibility<T>::operator() (size_t ik, size_t iw, size_t j, size_t jp, size_t m1,  size_t m2) const
{
	return *(this->read_phs_grid_ptr_block(ik,iw)+this->memory_layout_combined_notation_4pt_scalar_obj(j,jp,m1,m2));
}

template<typename T>
T GeneralizedSusceptibility<T>::operator()
	(size_t ik, size_t iw, size_t j, size_t jp, size_t l1, size_t l2, size_t l3, size_t l4) const
{
	return *(this->read_phs_grid_ptr_block(ik,iw)+this->memory_layout_4pt_scalar_obj(j,jp,l1,l2,l3,l4));
}

template<typename T>
T & GeneralizedSusceptibility<T>::operator()
	(size_t ik, size_t iw, size_t j, size_t jp, size_t l1, size_t l2, size_t l3, size_t l4)
{
	return *(this->write_phs_grid_ptr_block(ik,iw)+this->memory_layout_4pt_scalar_obj(j,jp,l1,l2,l3,l4));
}

template<typename T>
void GeneralizedSusceptibility<T>::AdiabaticUpscale_::set(
		bT maxScaling, size_t maxScaleResolutionSteps, size_t maxScaleSteps)
{
	this->reset();

	maxScaling_ = maxScaling;
	maxScaleResolutionSteps_ = maxScaleResolutionSteps;
	maxNScale_ = maxScaleSteps;
}

template<typename T>
void GeneralizedSusceptibility<T>::AdiabaticUpscale_::reset()
{
	scaling_ = bT(1);
	nscale_ = 0;
	softKVector_ = std::vector<bT>();
	softOrbitalVector_ = std::vector<T>();
}

template<typename T>
void GeneralizedSusceptibility<T>::AdiabaticUpscale_::get_soft_vector(
		std::vector<bT> & softKVector, std::vector<T> & softOrbitalVector) const
{
	softKVector = softKVector_;
	softOrbitalVector = softOrbitalVector_;
}

template<typename T>
void GeneralizedSusceptibility<T>::AdiabaticUpscale_::track_instability(
		std::vector<bT> const& softKVector, std::vector<T> const& softOrbitalVector )
{
	softKVector_ = softKVector;
	softOrbitalVector_ = softOrbitalVector;
	nscaleResolutionSteps_++;
	scaling_ = bT(1.0)+bT(maxScaling_-1.0)*bT(nscaleResolutionSteps_)/bT(maxScaleResolutionSteps_);
}

template<typename T>
typename GeneralizedSusceptibility<T>::bT
GeneralizedSusceptibility<T>::AdiabaticUpscale_::get_scaling() const
{
	return scaling_;
}

template<typename T>
void GeneralizedSusceptibility<T>::AdiabaticUpscale_::set_stable()
{
	if ( nscaleResolutionSteps_ > 0 )
	{
		nscale_++;
	}
	scaling_ = bT(1);
	nscaleResolutionSteps_ = 0;
}

template<typename T>
bool GeneralizedSusceptibility<T>::AdiabaticUpscale_::is_soft() const
{
	return nscaleResolutionSteps_ > 0;
}

template<typename T>
bool GeneralizedSusceptibility<T>::AdiabaticUpscale_::steps_ok() const
{
	return (nscaleResolutionSteps_ < maxScaleResolutionSteps_) && (nscale_ < maxNScale_);
}

template<typename T>
bool GeneralizedSusceptibility<T>::has_charge_part() const
{
	return channelsOffset_ == 0;
}

template<typename T>
size_t GeneralizedSusceptibility<T>::channels_offset() const
{
	return channelsOffset_;
}

} /* namespace gw_flex */
} /* namespace scallop */
