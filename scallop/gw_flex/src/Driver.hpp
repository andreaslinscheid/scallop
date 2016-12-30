/*	This file Driver.hpp is part of scallop.
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
 *  Created on: Dec 13, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/Driver.h"
#include "scallop/gw_flex/DysonEquation.h"
#include "scallop/input/KPath.h"
#include "scallop/parallel/IrregularGridDistribution.h"
#include "scallop/gw_flex/ProjectOutBandStructure.h"
#include "scallop/gw_flex/ManyBodyBandStructure.h"
#include "scallop/output/ObservableStatistics.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
Driver<T>::Driver(input::Configuration config) : config_(std::move(config))
{

}

template<typename T>
void Driver<T>::initialize( double temp, double & Ne)
{
	iter_ = 0;

	auto beta = auxillary::BasicFunctions::inverse_temperature( temp );

	elstr_.initialize_from_file( config_.get_grid() , config_.get_f_elstr() );
	if ( Ne > 0 )
	{
		elstr_.adjust_filling( Ne, beta );
	}
	else
	{
		output::TerminalOut msg(  auxillary::globals::VerbosityLvl::medium );
		Ne = elstr_.compute_N_electrons( beta );
		msg << "The number of electrons/U.C. is " << Ne;
	}

	auto mCutoff = config_.get_MCut();
	if ( mCutoff < 0)
	{
		output::TerminalOut msg(  auxillary::globals::VerbosityLvl::medium );
		auto bw = elstr_.get_band_width();
		mCutoff = std::abs(bw.first) > std::abs(bw.second) ? std::abs(bw.first) : std::abs(bw.second) ;
		mCutoff *= 1.5;
		msg << "Choosing a cutoff of " << mCutoff << " meV as 1.5 x max energy.";
	}

	size_t nM = 2*std::floor( mCutoff / ( 2.0*M_PI / beta ) - 0.5 );
	matsFreq_ = V(nM);
	for ( size_t n = 0 ; n < nM; ++n)
		matsFreq_[n] = auxillary::BasicFunctions::matzubara_frequency_of_index(n,nM,beta);

	if ( ! se_.is_init() )
	{
		V data;
		se_.initialize(
				matsFreq_.size(),
				config_.get_grid(),
				elstr_.get_nOrb(),
				/* init in time=*/ true,
				/* init in k spcae=*/ false,
				data );

		seL_ = se_;
	}
	else
	{
		//We possibly need to interpolate in frequency
		//and we need to update the number of electrons from a possible
		//non-zero self energy
	}
}

template<typename T>
void Driver<T>::converge()
{
	output::TerminalOut msg( auxillary::globals::VerbosityLvl::low );

	if ( not config_.get_f_sf().empty() )
	{
		msg << "Loading file " << config_.get_f_sf() << " with spin-fluctuation interaction parameters.";
		Isf_.init_file( config_.get_f_sf() );
	}

	if ( not config_.get_f_c().empty() )
	{
		msg << "Loading file " << config_.get_f_c() << " with charge-fluctuation interaction parameters.";
		Ic_.init_file( config_.get_f_c() );
	}

	for ( auto temp : config_.get_temp() )
	{
		for ( auto Ne : config_.get_nelec() )
		{
			output::TerminalOut msgDetailed( auxillary::globals::VerbosityLvl::high );
			std::string electrons;
			if ( Ne > 0 )
				electrons = std::string()+" with Ne=" + std::to_string(Ne)+" e/spin/U.C.  in the system";
			msg << "\nRunning at temperature T=" << temp << "K"+electrons+"\n";
			this->initialize(temp,Ne);

			msgDetailed << "Number of Matsubara frequencies is " << matsFreq_.size();
			while ( true )
			{
				iter_++;
				auto beta = auxillary::BasicFunctions::inverse_temperature( temp );
				msgDetailed << "\tStarting iteration " << iter_ << ":" ;
				this->prepare_iteration(temp);

				if ( config_.get_adj_mu() && iter_ > 1 )
				{
					msgDetailed << "\tAdjusting the chemical potential ...";
					this->shift_chemical_pot(temp,Ne);
				}
				G_.transform_itime_Mfreq_subtract(beta, ksGT_, ksGF_ );
				//Now G is in real space and time
				assert( G_.is_in_time_space() and (not G_.is_in_k_space()) );

				msgDetailed << "\tComputing the effective interactions ...";
				this->compute_eff_interactions( beta );

				msgDetailed << "\tConstruction self-energies ...";
				this->compute_self_energies();

				bool converged = this->check_convergence();

				this->mix_iterations();

				msgDetailed << "\tSolving Dyson's equation ...";
				se_.transform_itime_Mfreq( beta );
				se_.perform_space_fft( );
				this->solve_Dyson();


				assert( (! G_.is_in_time_space()) && (  G_.is_in_k_space()) );
				G_.perform_space_fft();
				this->report( beta );

				if ( converged )
				{
					msg << "Convergence achieved in " << iter_ << " iterations.";
					break;
				}

				if ( iter_ >= config_.get_maxIter() )
				{
					msg << "===============================================";
					msg << "| Warning: Convergence has not been achieved! |";
					msg << "===============================================";
					break;
				}
			}

			this->post_process();
		}
	}
}

template<typename T>
void Driver<T>::prepare_iteration( double temp)
{
	spinAdiabaticScale_.set(
			config_.get_maxDownscale(),
			config_.get_scale_res(),
			config_.get_scale_SE_iter());

	auto beta = auxillary::BasicFunctions::inverse_temperature( temp );
	//In the first iteration, we copy
	if ( (iter_ == 1) || ( config_.get_adj_mu() ) )
	{
		ksGF_.set_from_KS_bandstructure(
				/* in time space */ false,
				matsFreq_.size(),
				beta,
				elstr_ );
		ksGF_.perform_space_fft();

		ksGT_.set_from_KS_bandstructure(
				/* in time space */ true,
				matsFreq_.size(),
				beta,
				elstr_ );
		ksGT_.perform_space_fft();
	}

	//We start each iteration in space and frequency
	if (iter_ == 1)
	{
		G_ = ksGF_;
	}

	assert( (! G_.is_in_time_space())    && (! G_.is_in_k_space()) );
	assert( (! ksGF_.is_in_time_space()) && (! ksGF_.is_in_k_space()) );
	assert( (  ksGT_.is_in_time_space()) && (! ksGT_.is_in_k_space()) );
}

template<typename T>
void Driver<T>::shift_chemical_pot( double temp, double Ne  )
{
	std::pair<bool,bT> result;
	auto beta = auxillary::BasicFunctions::inverse_temperature( temp );
	chmpot_.determine_new_chemPot( Ne, beta, elstr_, ksGF_, G_, result);
	if ( result.first )
	{
		bT old_mu = elstr_.get_chem_pot();
		elstr_.set_chem_pot( result.second );
		ksGF_.set_chem_pot( result.second );
		G_.set_chem_pot( result.second );

		output::TerminalOut msg( auxillary::globals::VerbosityLvl::medium );
		msg << "\tShifting the chemical potential by " << result.second - old_mu << " meV";
	}
}

template<typename T>
void Driver<T>::compute_eff_interactions( bT beta )
{
	suscSF_.set_uninitialized();
	suscSF_.set_time_space( /*in time space =*/true );
	suscSF_.set_k_space( /*in time space =*/false );

	suscC_.set_uninitialized();
	suscC_.set_time_space( /*in time space =*/true );
	suscC_.set_k_space( /*in time space =*/false );

	if ( not Isf_.empty() )
	{
		suscSF_.compute_from_gf( G_ );
		if ( not Ic_.empty() )
		{
			suscC_.copy_charge_part( suscSF_ );
			suscC_.transform_itime_Mfreq( beta );
			suscC_.perform_space_fft( );
			suscC_.charge_RPA_enhancement( Ic_ );
			suscC_.transform_itime_Mfreq( beta );
			suscC_.perform_space_fft( );
		}
		suscSF_.transform_itime_Mfreq( beta );
		suscSF_.perform_space_fft( );
		suscSF_.spin_RPA_enhancement( Isf_, spinAdiabaticScale_ );
		while ( spinAdiabaticScale_.is_soft() )
		{
			suscSF_.set_uninitialized();
			suscSF_.set_time_space( /*in time space =*/true );
			suscSF_.set_k_space( /*in time space =*/false );
			suscSF_.compute_from_gf( G_ );
			suscSF_.transform_itime_Mfreq( beta );
			suscSF_.perform_space_fft( );
			suscSF_.spin_RPA_enhancement( Isf_, spinAdiabaticScale_ );

			if ( spinAdiabaticScale_.is_soft() and not spinAdiabaticScale_.steps_ok() )
			{
				//TODO magnetic, induce magnetism
				output::TerminalOut msg( auxillary::globals::VerbosityLvl::high );
				msg << "We are beyond the magnetic transition";
				std::abort();
			}
		}

		suscSF_.transform_itime_Mfreq( beta );
		suscSF_.perform_space_fft( );
	}
	else
	{
		if ( not Ic_.empty() )
		{
			suscC_.compute_from_gf( G_ );
			suscC_.transform_itime_Mfreq( beta );
			suscC_.perform_space_fft( );
			suscC_.charge_RPA_enhancement( Ic_ );
			suscC_.transform_itime_Mfreq( beta );
			suscC_.perform_space_fft( );
		}
	}

	assert( (suscC_.is_in_time_space() and (not suscC_.is_in_k_space())) );
	assert( (suscSF_.is_in_time_space() and (not suscSF_.is_in_k_space())) );
}

template<typename T>
void Driver<T>::compute_self_energies()
{
	se_.set_to_zero();
	se_.set_time_space( /*in time space =*/true );
	se_.set_k_space( /*in time space =*/false );
	if ( not Isf_.empty() )
		se_.add_electronic_selfenergy( G_ , suscSF_);
	if ( not Ic_.empty() )
		se_.add_electronic_selfenergy( G_ , suscC_);
}

template<typename T>
bool Driver<T>::check_convergence() const
{
	T cmplxZero = T(config_.get_cmp_zero(),config_.get_cmp_zero());
	return se_.compare_significant_digits(seL_,
			config_.get_acc_d(),
			cmplxZero);
}

template<typename T>
void Driver<T>::mix_iterations()
{
	//Make sure the Self-energies are in the same space
	assert( !(se_.is_in_k_space() xor seL_.is_in_k_space()) );
	assert( !(se_.is_in_time_space() xor seL_.is_in_time_space()) );

	auto mix = config_.get_mix();

	size_t nG = se_.is_in_k_space() ? se_.get_spaceGrid_proc().get_num_k_grid():se_.get_spaceGrid_proc().get_num_R_grid();
	for ( size_t ig = 0; ig < nG; ++ig )
		for ( size_t iw = 0; iw < se_.get_num_time(); ++iw )
		{
			auto ptrNew = se_.write_phs_grid_ptr_block(ig,iw);
			auto ptrOld = seL_.write_phs_grid_ptr_block(ig,iw);
			for ( size_t ib = 0 ; ib < se_.get_data_block_size(); ++ib)
			{
				ptrNew[ib] = ptrOld[ib]*bT(mix)+ptrNew[ib]*bT(1.0-mix);
				ptrOld[ib] = ptrNew[ib];
			}
		}
}

template<typename T>
void Driver<T>::solve_Dyson( )
{
	DysonEquation d;
	d.solve_by_inversion( G_, matsFreq_, elstr_, se_ );
}

template<typename T>
void Driver<T>::report( bT beta ) const
{
	output::TerminalOut msg( auxillary::globals::VerbosityLvl::medium );
	output::ObservableStatistics<T> o;
	o.print_statistics( msg, se_, elstr_ , G_, beta );
}

template<typename T>
void Driver<T>::post_process()
{
	output::TerminalOut msg( auxillary::globals::VerbosityLvl::medium );

	input::KPath<bT> path;
	if ( not config_.get_f_kpath().empty() )
	{
		path.read_kpath_file( config_.get_f_kpath() );
		if ( path.get_dim() != se_.get_spaceGrid_proc().get_dim() )
			error_handling::Error( std::string("Space dimension of the kpath in file '")+config_.get_f_kpath()
					+"' and the regular grid do not agree.");
	}

	if ( not config_.get_f_spec().empty() )
	{
		msg << "\tComputing the single particle spectral function ... ";
		ProjectOutBandStructure<T> bandProj( elstr_ );

		decltype(matsFreq_) mappedMatsFreq;
		typename auxillary::TemplateTypedefs<std::complex<float> >::scallop_vector mappedData;
		bandProj.project_out_KS_bands_half_freq(
				se_,
				matsFreq_,
				mappedMatsFreq, mappedData);

		bT minf, maxf, eta;
		size_t npts = 1000;
		if ( config_.get_r_omega().empty() )
		{
			auto bndw = elstr_.get_band_width();
			bT bw = std::abs(bndw.first - bndw.second);
			minf = bndw.first-bw*0.1;
			maxf = bndw.second+bw*0.1;
			eta = bw*1e-5;
		}
		else
		{
			assert( config_.get_r_omega().size() == 4 );
			minf = config_.get_r_omega()[0];
			maxf = config_.get_r_omega()[1];
			npts = static_cast<size_t>(std::floor(config_.get_r_omega()[2]+0.5));
			eta = config_.get_r_omega()[3];
		}

		decltype(matsFreq_) omega(npts);
		for (size_t io = 0 ;io  <npts ; ++io )
			omega[io] = T(minf + (bT(io)+0.5)/bT(npts)*(maxf-minf),0.5);

		parallel::IrregularGridDistribution<T> irrPath;
		irrPath.distribute_pts(
				/*in k space=*/ true,
				path.get_k_path(),
				se_.get_spaceGrid_proc());

		typename ProjectOutBandStructure<T>::container_type result;
		ManyBodyBandStructure<T> mb;
		mb.compute_spectral_function( irrPath,
					omega,
					true,
					se_.get_spaceGrid_proc(),
					mappedData,
					mappedMatsFreq,
					se_.get_data_block_size(),
					bandProj,
					result);

		path.band_structure_gnuplot( result,
				omega,
				irrPath,
				config_.get_f_spec() ,
				"${\\sum}_{n} A_{\\boldsymbol{k},n}$ [meV]");

		msg << "\tdone. Data and plotscript at "<<
				config_.get_f_spec() << ".dat, and "<<config_.get_f_spec()<<".gp, respectively";
	}
}

} /* namespace gw_flex */
} /* namespace scallop */
