/*	This file GreensFunctionOrbital_test.hpp is part of scallop.
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
 *  Created on: Nov 2, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/test/GreensFunctionOrbital_test.h"
#include "scallop/gw_flex/UnitaryWannierKSBands.h"
#include "scallop/output/TerminalOut.h"
#include "scallop/parallel/GridDistribution.h"
#include <cmath>
#include <fstream>

namespace scallop
{
namespace gw_flex
{
namespace test
{

template<typename T>
void GreensFunctionOrbital_test<T>::test_all()
{
	create_test_gfs();

	transform_reciprocal_to_realspace();

	test_full_loop_time_space_back();
}

template<typename T>
void GreensFunctionOrbital_test<T>::create_test_gfs()
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	output::TerminalOut msg;
	msg << "Initializing Fourier transform in k test of a GF of a simple cosine band model";

	//Here we define a single, 2D, cosine kx plus cosine ky band without orbital character
	const std::vector<size_t> kgrid = { 13 , 11 };
	const size_t nMKS = 128;
	const size_t nM = 256;
	const bT kb = 0.86173324 ; // meV / K
	const bT temp = 0.1; // K
	const bT beta = 1.0 / (kb * temp);
	const bT bandwidth = 20; // meV
	const bT mu = bandwidth*0.9;

	auto cos_bnd = [&] (size_t ikx, size_t iky, size_t sa)
		{
			return (sa < 2 ? 1.0:-1.0 )*(0.5*bandwidth*(std::cos( (2*M_PI*ikx)/kgrid[0] )
				+std::cos( (2*M_PI*iky)/kgrid[1] ))-mu);
		};

	parallel::GridDistribution<T> gd;
	gd.distribute_grid(kgrid);

	std::vector<bT> energies( gd.get_num_k_grid()*4 );
	for (size_t ik = 0 ; ik < gd.get_num_k_grid() ; ++ik)
	{
		auto tuple = gd.k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple.front();
		size_t iky = tuple.back();

		for (size_t sa = 0 ; sa < 4 ; ++sa)
			energies[ik*4+sa] = cos_bnd(ikx,iky,sa);
	}

	UnitaryWannierKSBands<T> unitary;
	unitary.initialize_identity(kgrid,1);

	//fill with the GF with data for the single band
	auto a = unitary;
	KS_gf_singleBand_.set_in_frequency_space( std::move(a), energies, nMKS, beta);
	T diff = T(0);
	for (size_t ik = 0 ; ik < KS_gf_singleBand_.get_spaceGrid_proc().get_num_k_grid() ; ++ik)
	{
		auto tuple = gd.k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple.front();
		size_t iky = tuple.back();

		for (size_t iw = 0 ; iw < nMKS; ++iw )
		{
			for (size_t sa = 0 ; sa < 4 ; ++sa)
			{
				int frequencyIndex =
						(iw < nMKS/2 ? static_cast<int>(iw) : static_cast<int>(iw)-static_cast<int>(nMKS) );

				bT energy = cos_bnd(ikx,iky,sa);
				T analytic = 1.0 / ( T(0,M_PI / beta * ( 2*frequencyIndex+1 ) ) - energy );

				diff += std::abs(std::real(KS_gf_singleBand_(ik,iw,sa,sa))-std::real(analytic));
				diff += std::abs(std::imag(KS_gf_singleBand_(ik,iw,sa,sa))-std::imag(analytic));
			}
		}
	}

	mpi.sum( diff );
	msg << "Difference between known form of the KS GF in frequency and the actual return of the object:" << diff;
	mpi.barrier();
	assert( (diff.real() < 0.0000001 ) && (diff.imag() < 0.0000001 ) );

	//fill with the GF with data for the single band
	KS_gf_singleBand_.set_in_time_space( std::move(unitary), std::move(energies), nMKS, beta);
	auto FermiFunc = [&] ( bT e) { return 1.0 / ( std::exp( beta * e ) + 1.0 );};
	diff = T(0);
	for (size_t ik = 0 ; ik < KS_gf_singleBand_.get_spaceGrid_proc().get_num_k_grid() ; ++ik)
	{
		auto tuple = gd.k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple.front();
		size_t iky = tuple.back();
		for (size_t iw = 0 ; iw < nMKS; ++iw )
		{
			for (size_t sa = 0 ; sa < 4 ; ++sa)
			{
				bT energy = cos_bnd(ikx,iky,sa);
				bT taui = (beta*iw) / nMKS;
				bT analytic = -1.0*FermiFunc(-energy)*std::exp(-taui*energy);

				diff += std::abs(std::real(KS_gf_singleBand_(ik,iw,sa,sa))-std::real(analytic));
				diff += std::abs(std::imag(KS_gf_singleBand_(ik,iw,sa,sa))-std::imag(analytic));
			}
		}
	}
	mpi.sum( diff );
	msg << "Difference between known form of the KS GF in time and the actual return of the object:" << diff;
	mpi.barrier();
	assert( (diff.real() < 0.0000001) && (diff.imag() < 0.0000001) );

	//initialize empty greens function for this grid
	typename auxillary::TemplateTypedefs<T>::scallop_vector gfd;
	gf_singleBand_.initialize(nM,kgrid,1,false,true,gfd);
}

template<typename T>
void GreensFunctionOrbital_test<T>::transform_reciprocal_to_realspace()
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	output::TerminalOut msg;
	msg << "Testing the reciprocal to direct lattice FFT of a cosine band.";

	std::vector<size_t> kgrid = gf_singleBand_.get_spaceGrid_proc().get_grid();
	const size_t nM = gf_singleBand_.get_num_time();
	const bT bandwidth = 20; // meV
	const bT mu = bandwidth*0.9;

	//For testing purposes, we set the GF data to t/4*(cos(kx)+cos(ky))
	//	since we know the Fourier transform analytically. Of cause this
	//	is not a real Green's function.
	for (size_t ik = 0 ; ik < KS_gf_singleBand_.get_spaceGrid_proc().get_num_k_grid() ; ++ik)
	{
		auto tuple = KS_gf_singleBand_.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple.front();
		size_t iky = tuple.back();
		for (size_t iw = 0 ; iw < nM ; ++iw)
			for (size_t sa = 0 ; sa < 4 ; ++sa)
				gf_singleBand_(ik,iw,sa,sa) = (sa < 2 ? 1.0:-1.0 )*
									0.5*bandwidth*(std::cos( (2*M_PI*ikx)/kgrid[0] )
										+std::cos( (2*M_PI*iky)/kgrid[1] ))-mu;
	}

	gf_singleBand_.perform_space_fft();

	T diff = T(0);
	for (size_t iR = 0 ; iR < KS_gf_singleBand_.get_spaceGrid_proc().get_num_R_grid() ; ++iR)
	{
		auto tuple = KS_gf_singleBand_.get_spaceGrid_proc().R_conseq_local_to_xyz_total( iR );
		size_t iRx = tuple.front();
		size_t iRy = tuple.back();
		for (size_t iw = 0 ; iw < nM ; ++iw)
		{
			T analytic = T(0);
			if ((iRx==iRy)&&(iRx==0))
				analytic = -mu;

			if ( ((iRx == 0)&&((iRy == 1)||(iRy == kgrid[1]-1)))
				 || ((iRy == 0)&&((iRx == 1)||(iRx == kgrid[0]-1))) )
				analytic = 0.25*bandwidth;

			diff += std::abs(std::real(gf_singleBand_(iR,iw,0,0))-std::real(analytic));
			diff += std::abs(std::imag(gf_singleBand_(iR,iw,0,0))-std::imag(analytic));
		}
	}
	mpi.sum( diff );
	msg << "The difference between analytic and numeric FFT from k to R is " << diff;
	mpi.barrier();
	assert( (diff.real() < 0.0000001) && (diff.imag() < 0.0000001) );

	gf_singleBand_.perform_space_fft();
	diff = T(0);
	for (size_t ik = 0 ; ik < KS_gf_singleBand_.get_spaceGrid_proc().get_num_k_grid() ; ++ik)
	{
		auto tuple = KS_gf_singleBand_.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple.front();
		size_t iky = tuple.back();
		for (size_t iw = 0 ; iw < nM ; ++iw)
		{
			T analytic = 0.5*bandwidth*(std::cos( (2*M_PI*ikx)/kgrid[0] )
				+std::cos( (2*M_PI*iky)/kgrid[1] ))-mu;

			diff += std::abs(std::real(gf_singleBand_(ik,iw,0,0))-std::real(analytic));
			diff += T(0,std::abs(std::imag(gf_singleBand_(ik,iw,0,0))-std::imag(analytic)));
		}
	}
	mpi.sum( diff );
	msg << "The difference between analytic and numeric FFT from R to k is " << diff;
	mpi.barrier();
	assert( (diff.real() < 0.0000001) && (diff.imag() < 0.0000001) );
}

template<typename T>
void GreensFunctionOrbital_test<T>::test_full_loop_time_space_back()
{
	output::TerminalOut msg;
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	msg << "Testing the time / frequency Fourier transform for a"
			"16x11 space grid for a single orbital free Green's function with a cosine dispersion at 500K with 2048 Matsubara points.";
	const std::vector<size_t> dimGrid = {16,11};
	parallel::GridDistribution<T> gd;
	gd.distribute_grid(dimGrid);

	const size_t nM = 2048;
	const bT kb = 0.86173324 ; // meV / K
	const bT temp = 500; // K
	const bT beta = 1.0 / (kb * temp);
	const bT bandwidth = 20; // meV
	const bT mu = bandwidth*0.9;

	const size_t orbitalDim = 1;
	const bool initalizeAsTime = false;
	const bool initalizeAsRecipr = true;
	typename auxillary::TemplateTypedefs<T>::scallop_vector gfd;
	gf_singleBand_.initialize(nM,dimGrid,orbitalDim,initalizeAsTime,initalizeAsRecipr,gfd);

	for (size_t ik = 0 ; ik < gf_singleBand_.get_spaceGrid_proc().get_num_k_grid() ; ++ik)
	{
		auto tuple = gf_singleBand_.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple.front();
		size_t iky = tuple.back();

		for (size_t i = 0 ; i < nM ; ++i)
			for (size_t as = 0 ; as < 4 ; ++as)
			{

				T energy = 0.5*bandwidth*(std::cos( (2*M_PI*ikx)/dimGrid[0] )
					+std::cos( (2*M_PI*iky)/dimGrid[1] ))-mu;

				int frequencyIndex =
						(i < nM/2 ? static_cast<int>(i) : static_cast<int>(i)-static_cast<int>(nM) );

				gf_singleBand_(ik,i,as,as) =
						1.0 / ( T(0,M_PI / beta * ( 2*frequencyIndex +1  ) ) - (as<2?1.0:-1.0)*energy );
			}
	}

	gf_singleBand_.transform_itime_Mfreq( beta );

	//compare the difference with the analytic formula
	auto fermiFunc = [] (bT energy, bT beta){ return 1.0 / ( std::exp( - beta * energy ) + 1.0 );};
	T diff = T(0);

	for (size_t ik = 0 ; ik < gf_singleBand_.get_spaceGrid_proc().get_num_k_grid() ; ++ik)
	{
		auto tuple = gf_singleBand_.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple.front();
		size_t iky = tuple.back();

		for (size_t i = 0 ; i < nM ; ++i)
			for (size_t as = 0 ; as < 4 ; ++as)
			{
				bT taui = (beta*i) / nM;
				bT energy = 0.5*bandwidth*(std::cos( (2*M_PI*ikx)/dimGrid[0] )
					+std::cos( (2*M_PI*iky)/dimGrid[1] ))-mu;
				T gAnalytic = -1.0*fermiFunc(energy,beta)*std::exp(-taui*energy);

				diff += std::abs( std::real(gf_singleBand_(ik,i,as,as))-std::real(gAnalytic))*beta/nM;
				diff += T(0,std::abs( std::imag(gf_singleBand_(ik,i,as,as))-std::imag(gAnalytic)))*(beta/nM);
			}
	}
	mpi.sum( diff );
	msg << "Difference between analytic and numeric Fourier transform from frequency to time of a free GF:"
			<< diff;
	mpi.barrier();
	assert( (diff.real() < 0.01) && (diff.imag() < 0.00000001) );


	msg << "Testing the time / frequency and space Fourier transform for a"
			"7x11x8 space grid for a single orbital free Green's function with a cosine dispersion at 500K with 2048 Matsubara points.";
	gd.distribute_grid( {7,11,8} );

	gf_singleBand_.initialize(nM,dimGrid,orbitalDim,initalizeAsTime,initalizeAsRecipr,gfd);

	for (size_t ik = 0 ; ik < gf_singleBand_.get_spaceGrid_proc().get_num_k_grid() ; ++ik)
	{
		auto tuple = gf_singleBand_.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple[0];
		size_t iky = tuple[1];
//		size_t ikz = tuple[2];//not used at the moment

		for (size_t i = 0 ; i < nM ; ++i)
			for (size_t as = 0 ; as < 4 ; ++as)
			{

				T energy = 0.5*bandwidth*(std::cos( (2*M_PI*ikx)/dimGrid[0] )
					+std::cos( (2*M_PI*iky)/dimGrid[1] ))-mu;

				int frequencyIndex =
						(i < nM/2 ? static_cast<int>(i) : static_cast<int>(i)-static_cast<int>(nM) );

				gf_singleBand_(ik,i,as,as) =
						1.0 / ( T(0,M_PI / beta * ( 2*frequencyIndex +1  ) ) - (as<2?1.0:-1.0)*energy );
			}
	}

	msg << "Transforming to time ...";
	gf_singleBand_.transform_itime_Mfreq( beta );

	msg << "Transforming to R ...";
	gf_singleBand_.perform_space_fft( );

	msg << "Transforming back to k ...";
	gf_singleBand_.perform_space_fft( );

	msg << "Cross-checking the difference to the analytic function";
	for (size_t ik = 0 ; ik < gf_singleBand_.get_spaceGrid_proc().get_num_k_grid() ; ++ik)
	{
		auto tuple = gf_singleBand_.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple[0];
		size_t iky = tuple[1];

		for (size_t i = 0 ; i < nM ; ++i)
			for (size_t as = 0 ; as < 4 ; ++as)
			{
				bT taui = (beta*i) / nM;
				bT energy = 0.5*bandwidth*(std::cos( (2*M_PI*ikx)/dimGrid[0] )
					+std::cos( (2*M_PI*iky)/dimGrid[1] ))-mu;
				T gAnalytic = -1.0*fermiFunc(energy,beta)*std::exp(-taui*energy);

				diff += std::abs( std::real(gf_singleBand_(ik,i,as,as))-std::real(gAnalytic))*beta/nM;
				diff += T(0,std::abs( std::imag(gf_singleBand_(ik,i,as,as))-std::imag(gAnalytic)))*(beta/nM);
			}
	}
	mpi.sum( diff );
	msg << "Difference between analytic and numeric Fourier transform from frequency to time of a free GF:"
			<< diff;
	mpi.barrier();
	assert( (diff.real() < 0.1) && (diff.imag() < 0.00000001) );

	msg << "Transforming back to frequency ...";
	gf_singleBand_.transform_itime_Mfreq( beta );

	msg << "Cross-checking the difference to the analytic function";
	for (size_t ik = 0 ; ik < gf_singleBand_.get_spaceGrid_proc().get_num_k_grid() ; ++ik)
	{
		auto tuple = gf_singleBand_.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple[0];
		size_t iky = tuple[1];

		for (size_t i = 0 ; i < nM ; ++i)
			for (size_t as = 0 ; as < 4 ; ++as)
			{
				bT energy = 0.5*bandwidth*(std::cos( (2*M_PI*ikx)/dimGrid[0] )
					+std::cos( (2*M_PI*iky)/dimGrid[1] ))-mu;
				int frequencyIndex =
						(i < nM/2 ? static_cast<int>(i) : static_cast<int>(i)-static_cast<int>(nM) );

				T gAnalytic = 1.0 / ( T(0,M_PI / beta * ( 2*frequencyIndex+1 ) ) - energy );

				diff += std::abs( std::real(gf_singleBand_(ik,i,as,as))-std::real(gAnalytic))*beta/nM;
				diff += T(0,std::abs( std::imag(gf_singleBand_(ik,i,as,as))-std::imag(gAnalytic)))*(beta/nM);
			}
	}
	mpi.sum( diff );
	msg << "Difference between analytic and numeric Fourier transform from frequency to time of a free GF:"
			<< diff;
	mpi.barrier();
	assert( (diff.real() < 0.1) && (diff.imag() < 0.0001) );
}

template<typename T>
template<class bandstructure>
GreensFunctionOrbital<T> GreensFunctionOrbital_test<T>::construct_free_time_gf(
		bT temperature, std::vector<size_t> spaceGrid, size_t nTimeSteps,
		size_t nBnd,
		bandstructure const& bnd)
{
	const bT kb = 0.86173324 ; // meV / K
	const bT beta = 1.0 / (kb * temperature);

	typename auxillary::TemplateTypedefs<T>::scallop_vector data;

	GreensFunctionOrbital<T> gf;
	gf.initialize(nTimeSteps,std::move(spaceGrid),nBnd,/* in time */true,/* in k */true, data );

	auto fermiFunc = [] (bT energy, bT beta) {return 1.0 / ( std::exp( beta * energy ) + 1.0 ) ; };

	for (size_t ik = 0 ; ik < gf.get_spaceGrid_proc().get_num_k_grid() ; ++ik)
	{
		auto tuple = gf.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		for (size_t i = 0 ; i < nTimeSteps ; ++i)
			for ( size_t l1 = 0 ; l1 < nBnd ; ++l1 )
				for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
					for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
						for ( size_t l2 = 0 ; l2 < nBnd ; ++l2 )
							{
								bT energy = bnd(tuple,l1,a1,s1,l2,a1,s1);
								bT taui = beta*(i+0.5) / nTimeSteps;
								bT analytic = energy > 0 ?
										-1.0*fermiFunc(-energy,beta)*std::exp(-taui*energy) :
										-1.0*fermiFunc(energy,beta)*std::exp((beta-taui)*energy);
								gf(ik,i,l1,a1,s1,l2,a1,s1) = analytic;
							}
	}
	return gf;
}

template<typename T>
template<class bandstructure>
GreensFunctionOrbital<T> GreensFunctionOrbital_test<T>::construct_gf_bnd(
		bT temperature, std::vector<size_t> spaceGrid, size_t nTimeSteps,
		size_t nBnd,
		bandstructure const& bnd)
{
	const bT kb = 0.86173324 ; // meV / K
	const bT beta = 1.0 / (kb * temperature);

	typename auxillary::TemplateTypedefs<T>::scallop_vector data;

	GreensFunctionOrbital<T> gf;
	gf.initialize(nTimeSteps,std::move(spaceGrid),nBnd,/* in time */false,/* in k */true, data );

	for (size_t ik = 0 ; ik < gf.get_spaceGrid_proc().get_num_k_grid() ; ++ik)
	{
		auto tuple = gf.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );

		for (size_t i = 0 ; i < nTimeSteps ; ++i)
			for ( size_t l1 = 0 ; l1 < nBnd ; ++l1 )
				for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
					for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
						for ( size_t l2 = 0 ; l2 < nBnd ; ++l2 )
							for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
								for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
								{
									bT energy = bnd(tuple,l1,a1,s1,l2,a2,s2);

									int frequencyIndex =
											(i < nTimeSteps/2 ? static_cast<int>(i) : static_cast<int>(i)-static_cast<int>(nTimeSteps) );

									gf(ik,i,l1,a1,s1,l2,a2,s2) = (a1==a2? 1.0:0.0 )*( s1==s2 ? 1.0:0.0 )
									     * 1.0 / ( T(0,M_PI / beta * ( 2*frequencyIndex+1 ) ) - energy );
								}
	}

	return gf;
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
