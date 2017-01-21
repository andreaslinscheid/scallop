/*	This file GeneralizedSusceptibility_test.hpp is part of scallop.
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

#include "scallop/gw_flex/SpinSusceptibility.h"
#include "scallop/gw_flex/ChargeSusceptibility.h"
#include "scallop/output/TerminalOut.h"
#include "scallop/parallel/MPIModule.h"
#include "scallop/gw_flex/test/GreensFunctionOrbital_test.h"
#include "scallop/gw_flex/test/InteractionMatrix_test.hpp"
#include <complex>

namespace scallop
{
namespace gw_flex
{
namespace test
{

template<typename T>
class GeneralizedSusceptibility_test
{
public:
	typedef typename scallop::auxillary::TypeMapComplex<T>::type bT;

	void test_all();

private:

	SpinSusceptibility<T> sust_;

	ChargeSusceptibility<T> csust_;

	void test_gf_construction();

	void test_enhancement();
};

template<typename T>
void GeneralizedSusceptibility_test<T>::test_all()
{
	output::TerminalOut msg;
	msg << "Testing the generalized susceptibility.";

	this->test_gf_construction();

	this->test_enhancement();
}

template<typename T>
void GeneralizedSusceptibility_test<T>::test_gf_construction()
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	output::TerminalOut msg;

	GreensFunctionOrbital<T> gf;

	auto pauli_0 = [] (size_t i , size_t j) {return i == j ? T(1.0) 					: T(0.0);};
	auto pauli_x = [] (size_t i , size_t j) {return i != j ? T(1.0) 					: T(0.0);};
	auto pauli_y = [] (size_t i , size_t j) {return i != j ? (i<j?T(0,1.0)	:T(0,-1.0))	: T(0.0);};
	auto pauli_z = [] (size_t i , size_t j) {return i == j ? (i==0?T(1.0)	:T(-1.0)) 	: T(0.0);};

	if ( mpi.get_nproc() == 1 )
	{
		msg << "Testing for 1 space and time point if the construction works";

		msg << "\nTesting the 1 orbital case. GF is set to combinations of  tau_x,y,z  sigma_0,x,y,z.";
		auto data = typename auxillary::TemplateTypedefs<T>::scallop_vector( 16 );

		msg << "Testing G =  tau z  sigma_0 (the charge channel)";
		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
						data[((a1*2+s1)*2+a2)*2+s2] = pauli_z(a1,a2)*pauli_0(s1,s2);

		gf.initialize(1, {1,1} , 1 , true, false, data );

		sust_.compute_from_gf( gf );

		for ( size_t j = 0; j < 4 ; ++j)
			for ( size_t jp = 0; jp < 4 ; ++jp)
				if ( not ( std::abs(( (j==jp)? T(-1.0/4.0) : T(0) ) - sust_(0,0,j,jp,0,0)  ) < 0.00000001 ) )
					error_handling::Error("Test failed");

		msg << "Testing G =  tau z sigma_x (the spin x channel)";
		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
						data[((a1*2+s1)*2+a2)*2+s2] = pauli_z(a1,a2)*pauli_x(s1,s2);

		gf.initialize(1, {1,1} , 1 , true, false, data );

		sust_.compute_from_gf( gf );

		for ( size_t j = 0; j < 4 ; ++j)
			for ( size_t jp = 0; jp < 4 ; ++jp)
				if ( not ( std::abs( -( (j==jp)? j<2? T(1):T(-1) : T(0) )/4.0 - sust_(0,0,j,jp,0,0)  ) < 0.00000001 ) )
					error_handling::Error("Test failed");

		msg << "Testing G =  tau z  sigma_y (the spin y channel)";
		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
						data[((a1*2+s1)*2+a2)*2+s2] = pauli_z(a1,a2)*pauli_y(s1,s2);

		gf.initialize(1, {1,1} , 1 , true, false, data );

		sust_.compute_from_gf( gf );

		for ( size_t j = 0; j < 4 ; ++j)
			for ( size_t jp = 0; jp < 4 ; ++jp)
				if ( not ( std::abs( -( (j==jp)? j%2==0? T(1):T(-1) : T(0) )/4.0 - sust_(0,0,j,jp,0,0)  ) < 0.00000001 ) )
					error_handling::Error("Test failed");

		msg << "Testing G =  tau z  sigma_z (the spin z channel)";
		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
						data[((a1*2+s1)*2+a2)*2+s2] = pauli_z(a1,a2)*pauli_z(s1,s2);

		gf.initialize(1, {1,1} , 1 , true, false, data );

		sust_.compute_from_gf( gf );

		for ( size_t j = 0; j < 4 ; ++j)
			for ( size_t jp = 0; jp < 4 ; ++jp)
				if ( not ( std::abs( -( (j==jp)? j==1||j==2? T(-1):T(1) : T(0) )/4.0 - sust_(0,0,j,jp,0,0)  ) < 0.00000001 ) )
					error_handling::Error("Test failed");

		msg << "Testing G = i tau y i sigma_y (the SC singlet channel)";
		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
						data[((a1*2+s1)*2+a2)*2+s2] = -1.0*pauli_y(a1,a2)*pauli_y(s1,s2);

		gf.initialize(1, {1,1} , 1 , true, false, data );

		sust_.compute_from_gf( gf );

		for ( size_t j = 0; j < 4 ; ++j)
			for ( size_t jp = 0; jp < 4 ; ++jp)
				if ( not ( std::abs( -( (j==jp)? j==0? T(-1):T(1) : T(0) )/4.0 - sust_(0,0,j,jp,0,0)  ) < 0.00000001 ) )
					error_handling::Error("Test failed");

		msg << "Testing G = i tau y (- sigma_z) (the SC triplet tx channel)";
		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
						data[((a1*2+s1)*2+a2)*2+s2] = -T(0,1.0)*pauli_y(a1,a2)*pauli_z(s1,s2);

		gf.initialize(1, {1,1} , 1 , true, false, data );

		sust_.compute_from_gf( gf );

		for ( size_t j = 0; j < 4 ; ++j)
			for ( size_t jp = 0; jp < 4 ; ++jp)
				if ( not ( std::abs( -( (j==jp)? j==1? T(-1):T(1) : T(0) )/4.0 - sust_(0,0,j,jp,0,0)  ) < 0.00000001 ) )
					error_handling::Error("Test failed");

		size_t nOrb = 2;
		msg << "\nTesting the "<<nOrb<<" orbital case";
		data = typename auxillary::TemplateTypedefs<T>::scallop_vector( 16*nOrb*nOrb );

		auto gfBandsVal = [] (size_t l1, size_t l2) {
			return ( l1 == 0 ? T(1.0) : T(2.0) )*( l2 == 0 ? T(1.0) : T(3.0) );
		};

		msg << "Testing G = tau z sigma_0 (the charge channel)";
		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
						for ( size_t l1 = 0 ; l1 < nOrb ; ++l1 )
							for ( size_t l2 = 0 ; l2 < nOrb ; ++l2 )
								data[((((a1*2+s1)*nOrb+l1)*2+a2)*2+s2)*nOrb+l2] =
										pauli_z(a1,a2)*pauli_0(s1,s2)
										*gfBandsVal(l1,l2);

		gf.initialize(1, {1,1} , 2 , true, false, data );

		sust_.set_uninitialized();
		sust_.compute_from_gf( gf );

		auto expectBandsVal = [&] (size_t l1, size_t l2, size_t l3, size_t l4) {
			return gfBandsVal(l1,l3)*gfBandsVal(l4,l2);
		};

		for ( size_t j = 0; j < 4 ; ++j)
			for ( size_t jp = 0; jp < 4 ; ++jp)
				for ( size_t l1 = 0 ; l1 < nOrb ; ++l1 )
					for ( size_t l2 = 0 ; l2 < nOrb ; ++l2 )
						for ( size_t l3 = 0 ; l3 < nOrb ; ++l3 )
							for ( size_t l4 = 0 ; l4 < nOrb ; ++l4 )
							{
								auto v = expectBandsVal(l1,l2,l3,l4);
								if ( not ( std::abs( -( (j==jp)? v : T(0) )/4.0 - sust_(0,0,j,jp,l1*nOrb+l2,l3*nOrb+l4)  ) < 0.00000001 ) )
									error_handling::Error("Test failed");
							}
	}

	msg << "Testing for a grid of 4x4 with 2 time points and 1 band:";
	msg << "\tTesting G = tau z sigma_0 (the charge channel)";
	//We are testing the construction by making the second time value imaginary and
	//adding a prefactor of the consecutive grid index. Thus, the susceptibility must
	//be purely imaginary and proportional to the consequtive grid index squared

	typename auxillary::TemplateTypedefs<T>::scallop_vector data;
	gf.initialize(2, {4,4}, 1, true, false, data );

	for ( size_t iR = 0 ; iR < gf.get_spaceGrid_proc().get_num_R_grid(); ++iR)
	{
		auto tuple = gf.get_spaceGrid_proc().R_conseq_local_to_xyz_total( iR );
		auto grid = gf.get_spaceGrid_proc().get_grid();
		for ( size_t it = 0 ; it < gf.get_num_time(); ++it)
			for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
				for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
					for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
						for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
						{
							T pref = (tuple[0]+1)+tuple[1]*grid[0];
							pref = (it == 0? pref : T(0,1.0)*pref);
							gf(iR,it,0,a1,s1,0,a2,s2) = pref*pauli_z(a1,a2)*pauli_0(s1,s2);
						}
	}
	sust_.set_uninitialized();
	sust_.compute_from_gf( gf );

	for ( size_t iR = 0 ; iR < gf.get_spaceGrid_proc().get_num_R_grid(); ++iR)
	{
		auto tuple = gf.get_spaceGrid_proc().R_conseq_local_to_xyz_total( iR );
		auto grid = gf.get_spaceGrid_proc().get_grid();
		for ( size_t it = 0 ; it < gf.get_num_time(); ++it)
			for ( size_t j = 0; j < 4 ; ++j)
				for ( size_t jp = 0; jp < 4 ; ++jp)
				{
					T pref = tuple[0]+tuple[1]*grid[0]+1;
					auto v =  T(0,1.0)*pref*pref;
					if ( not (std::abs( -( (j==jp)? v : T(0) )/4.0 - sust_(iR,it,j,jp,0,0)  ) < 0.00000001) )
						error_handling::Error("Test failed");
				}
	}

	msg << "Testing to get the bare suscepibility of a cos kx + cos ky band model.";
	const bT bandwidth = 20; // meV
	const bT mu = 0; // Perfect nesting
	const bT temperature = 0.3; // Temperature hand tuned to yield -0.1 to good accuracy at Q=pi,pi
	const bT beta = auxillary::BasicFunctions::inverse_temperature(temperature) ;
	const size_t nM = 128;

	std::vector<size_t> spaceGrid = { 64, 64 };
	msg << "\tThe grid is " << spaceGrid[0] << "x"<< spaceGrid[1] << " with "<< nM << " time steps (Matsubara points).";

	auto cos_bnd = [&] (std::vector<size_t> const& tuple,
			size_t l1, size_t a1, size_t s1,
			size_t l2, size_t a2, size_t s2){
		return ( a1 == a2 ? a1 == 0 ? 1.0 : -1.0 : 0.0 ) *  ( s1 == s2 ? 1.0 : 0.0 )
				*( 0.5*bandwidth*(std::cos( (2*M_PI*tuple[0])/spaceGrid[0] )
							+std::cos( (2*M_PI*tuple[1])/spaceGrid[1] ))-mu);
	};

	msg << "\tConstructing green's function in k and time domain ...";
	GreensFunctionOrbital_test<T> gftest;
	gf = gftest.construct_free_time_gf(
			temperature, spaceGrid ,
			nM, 1,
			cos_bnd);

	msg << "\tFT to lattice space ...";
	gf.perform_space_fft();

	msg << "\tConstructing susceptibility ...";
	sust_.set_uninitialized();
	sust_.compute_from_gf( gf );

	msg << "\tFT to reciprocal space ...";
	sust_.perform_space_fft();

	msg << "\tFT to frequency space ...";
	sust_.transform_itime_Mfreq(beta);

	size_t mid = spaceGrid[1]*(spaceGrid[0]/2)+spaceGrid[1]/2;
	for (size_t ik = 0 ; ik < sust_.get_spaceGrid_proc().get_num_k_grid(); ++ik)
	{
		auto tuple = sust_.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ictotal = sust_.get_spaceGrid_proc().k_xyz_to_conseq( tuple );
		if ( ictotal == mid)
		{
			std::cout << "\tValue at Q=(pi,pi):" << sust_(ik,0,0,0,0,0) << std::endl;
			if ( not ( (std::abs( std::real(sust_(ik,0,0,0,0,0))+0.42994) < 0.00001)
					and (std::abs( std::imag(sust_(ik,0,0,0,0,0)) ) < 0.00001 ) ) )
				error_handling::Error("Test failed");
		}
		if ( ictotal == 0 )
		{
			std::cout << "\tValue at Q=(0,0):" << sust_(ik,0,0,0,0,0) << std::endl;
		}
	}
}

template<typename T>
void GeneralizedSusceptibility_test<T>::test_enhancement()
{
	output::TerminalOut msg;

	msg << "Testing the copy of the charge channel";
	csust_.copy_charge_part(sust_);

	auto spaceGrid = csust_.get_spaceGrid_proc().get_grid();
	size_t mid = spaceGrid[1]*(spaceGrid[0]/2)+spaceGrid[1]/2;
	for (size_t ik = 0 ; ik < csust_.get_spaceGrid_proc().get_num_k_grid(); ++ik)
	{
		auto tuple = csust_.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ictotal = csust_.get_spaceGrid_proc().k_xyz_to_conseq( tuple );
		if ( ictotal == mid)
		{
			std::cout << "\tValue at Q=(pi,pi):" << csust_(ik,0,0,0,0,0) << std::endl;
			assert( (std::abs( std::real(csust_(ik,0,0,0,0,0))+0.42994) < 0.00001)
					and (std::abs( std::imag(csust_(ik,0,0,0,0,0)) ) < 0.00001 ) );
		}
		if ( ictotal == 0 )
		{
			std::cout << "\tValue at Q=(0,0):" << csust_(ik,0,0,0,0,0) << std::endl;
		}
	}

	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	const std::string fnameSpin = "/tmp/test_input_s.dat";
	if ( mpi.ioproc() )
	{
		std::ofstream testFile( fnameSpin.c_str() );
		//Testing InteractionMatrix input from file:
		testFile << "1 4" << '\n'; //One orbital 4 channels
		testFile << "0\t0\t0\t0\t0\t0\t0.125\t0.0" << '\n';
		testFile << "1\t1\t0\t0\t0\t0\t0.25\t0.0" << '\n';
		testFile << "2\t2\t0\t0\t0\t0\t0.375\t0.0" << '\n';
		testFile << "3\t3\t0\t0\t0\t0\t0.5\t0.0" << '\n';
		//			^j ^jp l1 l2 l3 l4 Re(I) Im(I)
		testFile.close();
	}

	const std::string fnameCharge = "/tmp/test_input_c.dat";
	if ( mpi.ioproc() )
	{
		std::ofstream testFile( fnameCharge.c_str() );
		//Testing InteractionMatrix input from file:
		testFile << "1 1" << '\n'; //One orbital one channel
		testFile << "0\t0\t0\t0\t0\t0\t1.5\t0.0" << '\n';
		//			^j ^jp l1 l2 l3 l4 Re(I) Im(I)
		testFile.close();
	}
	mpi.barrier();

	InteractionMatrix<T> interactCharge;
	interactCharge.init_file( "/tmp/test_input_c.dat" );

	msg << "Testing the RPA charge enhancement";
	csust_.charge_RPA_enhancement( interactCharge );
	for (size_t ik = 0 ; ik < csust_.get_spaceGrid_proc().get_num_k_grid(); ++ik)
	{
		auto tuple = csust_.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ictotal = csust_.get_spaceGrid_proc().k_xyz_to_conseq( tuple );
		if ( ictotal == mid)
		{
			std::cout << "\tValue at Q=(pi,pi):" << csust_(ik,0,0,0,0,0) << std::endl;
//			assert( (std::abs( std::real(csust_(ik,0,0,0,0,0))-11.154) < 0.01)
//					and (std::abs( std::imag(csust_(ik,0,0,0,0,0)) ) < 0.001 ) );
		}
		if ( ictotal == 0 )
		{
			std::cout << "\tValue at Q=(0,0):" << csust_(ik,0,0,0,0,0) << std::endl;
		}
	}

	InteractionMatrix<T> interactSpin;
	interactSpin.init_file( "/tmp/test_input_s.dat" );

	msg << "Testing the RPA spin enhancement";
	sust_.spin_RPA_enhancement( interactSpin );
	for (size_t ik = 0 ; ik < sust_.get_spaceGrid_proc().get_num_k_grid(); ++ik)
	{
		auto tuple = sust_.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ictotal = sust_.get_spaceGrid_proc().k_xyz_to_conseq( tuple );
		if ( ictotal == mid)
		{
			std::cout << "\tValues at Q=(pi,pi):\n\tjp=\t1\t2\t3\t4";
			for (size_t j = 0 ; j < 4; ++j)
			{
				std::cout << "\n\tj=" << j;
				for (size_t jp = 0 ; jp < 4; ++jp)
				{
					std::cout << '\t' << sust_(ik,0,j,jp,0,0);
				}
			}
			std::cout << std::endl;
//			assert( (std::abs( std::real(sust_(ik,0,3,3,0,0))-10.5537) < 0.01)
//					and (std::abs( std::imag(sust_(ik,0,3,3,0,0)) ) < 0.001 ) );
		}
		if ( ictotal == 0 )
		{
			std::cout << "\tValue at Q=(0,0):\n\tj=\t1\t2\t3\t4";
			for (size_t j = 0 ; j < 4; ++j)
			{
				std::cout << "\n\tj=" << j;
				for (size_t jp = 0 ; jp < 4; ++jp)
				{
					std::cout << '\t' << sust_(ik,0,j,jp,0,0);
				}
			}
			std::cout << std::endl;
//			assert( (std::abs( std::real(sust_(ik,0,3,3,0,0))-0.265168) < 0.01)
//					and (std::abs( std::imag(sust_(ik,0,3,3,0,0)) ) < 0.001 ) );
		}
	}
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
