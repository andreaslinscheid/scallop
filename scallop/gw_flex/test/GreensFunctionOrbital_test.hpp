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
#include <cmath>
#include <fstream>

namespace scallop
{
namespace gw_flex
{
namespace test
{

template<typename T>
GreensFunctionOrbital_test<T>::GreensFunctionOrbital_test()
{
	std::cout << "Initializing Fourier transform in k test of a GF of a simple cosine band model" << std::endl;
	//Here we define a single, 2D, cosine kx plus cosine ky band without orbital character
	const std::vector<size_t> kgrid = { 32 , 32 };
	const size_t nMKS = 16;
	const size_t nM = 2;
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

	//initialize empty greens function for this grid
	std::vector<T> gfd(nM*16*kgrid[0]*kgrid[1], T(0) );
	gf_singleBand_.initialize(nM,kgrid,1,true,true,std::move(gfd));

	std::vector<bT> energies( kgrid[0]*kgrid[1]*4 );
	for (size_t ikx = 0 ; ikx < kgrid[0] ; ++ikx)
		for (size_t iky = 0 ; iky < kgrid[1] ; ++iky)
			for (size_t sa = 0 ; sa < 4 ; ++sa)
				energies[(ikx*kgrid[1]+iky)*4+sa] = cos_bnd(ikx,iky,sa);

	UnitaryWannierKSBands<T> unitary;
	unitary.initialize_identity(kgrid,1);

	//fill with the GF with data for the single band
	KS_gf_singleBand_.set_in_frequency_space( unitary,	energies, nMKS, beta);
	T diff = T(0);
	for (size_t ikx = 0 ; ikx < kgrid[0] ; ++ikx)
	{
		for (size_t iky = 0 ; iky < kgrid[1] ; ++iky)
		{
			for (size_t iw = 0 ; iw < nMKS; ++iw )
			{
				for (size_t sa = 0 ; sa < 4 ; ++sa)
				{
					int frequencyIndex =
							(iw < nMKS/2 ? static_cast<int>(iw) : static_cast<int>(iw)-static_cast<int>(nMKS) );

					bT energy = cos_bnd(ikx,iky,sa);
					T analytic = 1.0 / ( T(0,M_PI / beta * ( 2*frequencyIndex+1 ) ) - energy );

					diff += std::abs(std::real(KS_gf_singleBand_(ikx*kgrid[1]+iky,iw,sa,sa))-std::real(analytic));
					diff += std::abs(std::imag(KS_gf_singleBand_(ikx*kgrid[1]+iky,iw,sa,sa))-std::imag(analytic));
				}
			}
		}
	}
	std::cout << "Difference between known form of the KS GF in frequency and the actual return of the object:" << diff << std::endl;

	//fill with the GF with data for the single band
	KS_gf_singleBand_.set_in_time_space( unitary,	std::move(energies), nMKS, beta);
	auto FermiFunc = [&] ( bT e) { return 1.0 / ( std::exp( beta * e ) + 1.0 );};
	diff = T(0);
	for (size_t ikx = 0 ; ikx < kgrid[0] ; ++ikx)
	{
		for (size_t iky = 0 ; iky < kgrid[1] ; ++iky)
		{
			for (size_t iw = 0 ; iw < nMKS; ++iw )
			{
				for (size_t sa = 0 ; sa < 4 ; ++sa)
				{
					bT energy = cos_bnd(ikx,iky,sa);
					T analytic = -FermiFunc(-energy)*std::exp(-beta*energy);

					diff += std::abs(std::real(KS_gf_singleBand_(ikx*kgrid[1]+iky,iw,sa,sa))-std::real(analytic));
					diff += std::abs(std::imag(KS_gf_singleBand_(ikx*kgrid[1]+iky,iw,sa,sa))-std::imag(analytic));
				}
			}
		}
	}
	std::cout << "Difference between known form of the KS GF in time and the actual return of the object:" << diff << std::endl;
}

template<typename T>
void GreensFunctionOrbital_test<T>::transform_reciprocal_to_realspace()
{
	std::cout << "Testing the reciprocal to direct lattice FFT of a cosine band." << std::endl;

	const std::vector<size_t> kgrid = gf_singleBand_.get_spaceGrid_proc();
	const size_t nM = gf_singleBand_.get_num_time();
	const bT bandwidth = 20; // meV
	const bT mu = bandwidth*0.9;

	//For testing purposes, we set the GF data to t/4*(cos(kx)+cos(ky))
	//	since we know the Fourier transform analytically. Of cause this
	//	is not a real Green's function.
	for (size_t ikx = 0 ; ikx < kgrid[0] ; ++ikx)
		for (size_t iky = 0 ; iky < kgrid[1] ; ++iky)
			for (size_t iw = 0 ; iw < nM ; ++iw)
				for (size_t sa = 0 ; sa < 4 ; ++sa)
					gf_singleBand_(ikx*kgrid[1]+iky,iw,sa,sa) = (sa < 2 ? 1.0:-1.0 )*
												0.5*bandwidth*(std::cos( (2*M_PI*ikx)/kgrid[0] )
													+std::cos( (2*M_PI*iky)/kgrid[1] ))-mu;

	gf_singleBand_.perform_k_to_R_fft();

	T diff = T(0);
	for (size_t iRx = 0 ; iRx < kgrid[0] ; ++iRx)
	{
		for (size_t iRy = 0 ; iRy < kgrid[1] ; ++iRy)
		{
			for (size_t iw = 0 ; iw < nM ; ++iw)
			{
				T analytic = T(0);
				if ((iRx==iRy)&&(iRx==0))
					analytic = -mu;

				if ( ((iRx == 0)&&((iRy == 1)||(iRy == kgrid[1]-1)))
					 || ((iRy == 0)&&((iRx == 1)||(iRx == kgrid[0]-1))) )
					analytic = 0.25*bandwidth;

				diff += std::abs(std::real(gf_singleBand_(iRx*kgrid[1]+iRy,iw,0,0))-std::real(analytic));
				diff += std::abs(std::imag(gf_singleBand_(iRx*kgrid[1]+iRy,iw,0,0))-std::imag(analytic));

			}
		}
	}
	std::cout << "The difference between analytic and numeric FFT from k to R is " << diff << std::endl;

	gf_singleBand_.perform_R_to_k_fft();
	diff = T(0);
	for (size_t ikx = 0 ; ikx < kgrid[0] ; ++ikx)
	{
		for (size_t iky = 0 ; iky < kgrid[1] ; ++iky)
		{
			for (size_t iw = 0 ; iw < nM ; ++iw)
			{
				T analytic = 0.5*bandwidth*(std::cos( (2*M_PI*ikx)/kgrid[0] )
					+std::cos( (2*M_PI*iky)/kgrid[1] ))-mu;

				diff += std::abs(std::real(gf_singleBand_(ikx*kgrid[1]+iky,iw,0,0))-std::real(analytic));
				diff += std::abs(std::imag(gf_singleBand_(ikx*kgrid[1]+iky,iw,0,0))-std::imag(analytic));
			}
		}
	}
	std::cout << "The difference between analytic and numeric FFT from R to k is " << diff << std::endl;
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
