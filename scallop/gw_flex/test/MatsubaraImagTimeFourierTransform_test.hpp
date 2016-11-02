/*	This file MatsubaraImagTimeFourierTransform_test.hpp is part of scallop.
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
 *  Created on: Oct 28, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/GreensFunctionOrbital.h"
#include "scallop/auxillary/TypeMapComplex.h"
#include <vector>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <fstream>

namespace scallop
{
namespace gw_flex
{
namespace test
{

template<typename T>
void
MatsubaraImagTimeFourierTransform_test<T>::test_free_particle_greensfunction()
{
	typedef typename scallop::auxillary::TypeMapComplex<T>::type bT;

	//Minimal test for 1 grid point in kx and ky and one orbital
	const size_t nM = 10000;
	const bT kb = 0.86173324 ; // meV / K
	const bT temp = 0.1; // K
	const bT beta = 1.0 / (kb * temp);
	const bT energy = 10; // meV

	std::vector<T> gfd(nM*16, T(0) );
	for (size_t i = 0 ; i < nM ; ++i)
		for (size_t a = 0 ; a < 2 ; ++a)
			for (size_t s = 0 ; s < 2 ; ++s)
			{
				int frequencyIndex =
						(i < nM/2 ? static_cast<int>(i) : static_cast<int>(i)-static_cast<int>(nM) );

				gfd[(((i*2+a)*2+s)*2+a)*2+s] =
						1.0 / ( T(0,M_PI / beta * ( 2*frequencyIndex +1  ) ) - (a==0?1.0:-1.0)*energy );
			}

	scallop::gw_flex::GreensFunctionOrbital< std::complex<double> > gf;
	const size_t orbitalDim = 1;
	const std::vector<size_t> dimGrid = {1,1};
	const bool initalizeAsFreq = true;
	const bool initalizeAsRecipr = true;
	gf.initialize(nM,dimGrid,orbitalDim,initalizeAsFreq,initalizeAsRecipr,std::move(gfd));

	T diff = T(0);
	for (size_t i = 0 ; i < nM ; ++i)
	{
		int frequencyIndex =
				(i < nM/2 ? static_cast<int>(i) : static_cast<int>(i)-static_cast<int>(nM) );

		T ag = 1.0 / ( T(0,M_PI / beta * ( 2*frequencyIndex+1 ) ) - energy );
		diff += std::abs( std::real(gf(0,i,0,0,0,0,0,0))-std::real(ag));
		diff += T(0,std::abs( std::imag(gf(0,i,0,0,0,0,0,0))-std::imag(ag)));
	}
	std::cout << "Difference between the free GF and its initialized data:"
			<< diff << std::endl;

	gf.transform_Mfreq_to_itime( beta );

	//compare the difference with the analytic formula
	T fermiFunc =  1.0 / ( std::exp( - beta * energy ) + 1.0 );
	diff = T(0);
	for (size_t i = 0 ; i < nM ; ++i)
	{
		bT taui = (beta*i) / nM;
		T gAnalytic = -fermiFunc*std::exp(-taui*energy);
		diff += std::abs( std::real(gf(0,i,0,0,0,0,0,0))-std::real(gAnalytic))*beta/nM;
		diff += T(0,std::abs( std::imag(gf(0,i,0,0,0,0,0,0))-std::imag(gAnalytic)))*(beta/nM);
	}

	std::cout << "Difference between analytic and numeric Fourier transform from frequency to time of a free GF:"
			<< diff << std::endl;

	//Now, the other way around
	for (size_t i = 0 ; i < nM ; ++i)
		for (size_t a = 0 ; a < 2 ; ++a)
			for (size_t s = 0 ; s < 2 ; ++s)
			{
				T fermiFunc =  1.0 / ( std::exp( - beta * (a==0?1.0:-1.0)*energy ) + 1.0 );
				bT taui = (beta*i) / nM;
				gf(0,i,0,a,s,0,a,s) = -fermiFunc*std::exp(-taui*(a==0?1.0:-1.0)*energy);
			}
	gf.set_time_domain();
	gf.transform_itime_to_Mfreq( beta );
	diff = T(0);
std::ofstream file("/home/alinsch/codes/scallop/tests/bla2");
	for (size_t i = 0 ; i < nM ; ++i)
	{
		int frequencyIndex =
				(i < nM/2 ? static_cast<int>(i) : static_cast<int>(i)-static_cast<int>(nM) );

		T ag = 1.0 / ( T(0,M_PI / beta * ( 2*frequencyIndex+1 ) ) - energy );
		diff += std::abs( std::real(gf(0,i,0,0,0,0,0,0))-std::real(ag))*beta/nM;
		diff += T(0,std::abs( std::imag(gf(0,i,0,0,0,0,0,0))-std::imag(ag)))*(beta/nM);
file << i <<'\t'<< std::real(gf(0,i,0,0,0,0,0,0)) <<'\t'<< std::imag(gf(0,i,0,0,0,0,0,0))
			<<'\t'<< std::real(ag) <<'\t'<< std::imag(ag) <<'\n';
	}

	std::cout << "Difference between analytic and numeric Fourier transform from time to frequency of a free GF:"
			<< diff << std::endl;
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
