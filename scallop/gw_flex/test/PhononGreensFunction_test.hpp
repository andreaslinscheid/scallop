/*	This file PhononGreensFunction_test.hpp is part of scallop.
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
 *  Created on: Nov 14, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/PhononGreensFunction.h"
#include "scallop/output/TerminalOut.h"

namespace scallop
{
namespace gw_flex
{
namespace test
{

template<typename T>
class PhononGreensFunction_test
{
public:
	void test_all();

private:
	void test_with_Einstein_mode();

};

template<typename T>
void PhononGreensFunction_test<T>::test_all()
{
	test_with_Einstein_mode();
}

template<typename T>
void PhononGreensFunction_test<T>::test_with_Einstein_mode()
{
	typedef typename PhononGreensFunction<T>::bT bT;

	output::TerminalOut msg;
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	size_t nM = 2048;
	msg << "Testing initialization of a phonon propagator with a single Einstein mode on a 10x10 grid with "<<nM<<" Matsubara points";

	const bT modeEnergy = 10;
	const bT kb = 0.86173324 ; // meV / K
	const bT temperature = 1;
	const bT beta = 1.0 / (kb * temperature);

	PhononGreensFunction<T> pgf;
	typename auxillary::TemplateTypedefs<T>::scallop_vector data;
	pgf.initialize(nM, {10, 10}, 1, /*in time=*/ false, /*in k space=*/true, data );

	for ( size_t iq = 0 ; iq <pgf.get_spaceGrid_proc().get_num_k_grid(); ++iq )
	{
		for ( size_t iw = 0 ; iw < nM; ++iw )
		{
			int frequencyIndex =
					(iw < nM/2 ? static_cast<int>(iw) : static_cast<int>(iw)-static_cast<int>(nM) );

			pgf(iq,iw,0,0) = 2*modeEnergy / ( std::pow(2*M_PI*frequencyIndex/beta ,2) + std::pow(modeEnergy ,2) );
		}
	}

	msg << "Transforming to the time domain ...";
	pgf.transform_itime_Mfreq( beta );

	msg << "Comparing with the analytical result ...";
	T diff = T(0);
	for ( size_t iq = 0 ; iq <pgf.get_spaceGrid_proc().get_num_k_grid(); ++iq )
	{
		for ( size_t iw = 0 ; iw < nM; ++iw )
		{
			bT itau = beta * bT(iw+0.5) / bT(nM);
			T analytic = std::cosh( modeEnergy * ( 0.5*beta - itau ))/std::sinh(modeEnergy*0.5*beta);

			diff += std::abs(std::real(pgf(iq,iw,0,0))-std::real(analytic));
			diff += T(0,std::abs(std::imag(pgf(iq,iw,0,0))-std::imag(analytic)));
		}
	}
	diff *= (beta/nM);

	mpi.sum( diff );
	msg << "Difference between analytic and numeric Fourier transform from frequency to time:"
			<< diff;
	mpi.barrier();
	assert( (diff.real() < 0.001) && (diff.imag() < 0.0001) );

	msg << "Transforming back to the frequency domain ...";
	pgf.transform_itime_Mfreq( beta );
	msg << "Comparing with the analytical result ...";
	diff = T(0);
	for ( size_t iq = 0 ; iq <pgf.get_spaceGrid_proc().get_num_k_grid(); ++iq )
	{
		for ( size_t iw = 0 ; iw < nM; ++iw )
		{
			int frequencyIndex =
					(iw < nM/2 ? static_cast<int>(iw) : static_cast<int>(iw)-static_cast<int>(nM) );

			T analytic = 2*modeEnergy / ( std::pow(2*M_PI*frequencyIndex/beta ,2) + std::pow(modeEnergy ,2) );

			diff += std::abs(std::real(pgf(iq,iw,0,0))-std::real(analytic));
			diff += T(0,std::abs(std::imag(pgf(iq,iw,0,0))-std::imag(analytic)));
		}
	}
	diff *= (beta/nM);
	mpi.sum( diff );
	msg << "Difference between analytic and numeric Fourier transform from time to frequency:"
			<< diff;
	mpi.barrier();
	assert( (diff.real() < 0.0001) && (diff.imag() < 0.00001) );
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
