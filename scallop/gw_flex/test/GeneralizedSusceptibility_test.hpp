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

#include "scallop/gw_flex/GeneralizedSusceptibility.h"
#include "scallop/output/TerminalOut.h"
#include "scallop/parallel/MPIModule.h"
#include "scallop/gw_flex/test/GreensFunctionOrbital_test.h"
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

	GeneralizedSusceptibility<T> sust_;

	void test_gf_construction();
};

template<typename T>
void GeneralizedSusceptibility_test<T>::test_all()
{
	output::TerminalOut msg;
	msg << "Testing the generalized susceptibility.";

	this->test_gf_construction();
}

template<typename T>
void GeneralizedSusceptibility_test<T>::test_gf_construction()
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	output::TerminalOut msg;

	GreensFunctionOrbital<T> gf;

	if ( mpi.get_nproc() == 1 )
	{
		msg << "Testing for 1 space and time point if the construction works";

		msg << "\nTesting the 1 orbital case. GF is set to combinations of  tau_x,y,z  sigma_0,x,y,z.";
		auto data = typename auxillary::TemplateTypedefs<T>::scallop_vector( 16 );

		auto pauli_0 = [] (size_t i , size_t j) {return i == j ? T(1.0) 					: T(0.0);};
		auto pauli_x = [] (size_t i , size_t j) {return i != j ? T(1.0) 					: T(0.0);};
		auto pauli_y = [] (size_t i , size_t j) {return i != j ? (i<j?T(0,1.0)	:T(0,-1.0))	: T(0.0);};
		auto pauli_z = [] (size_t i , size_t j) {return i == j ? (i==0?T(1.0)	:T(-1.0)) 	: T(0.0);};

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
				assert( std::abs(( (j==jp)? T(-1.0/4.0) : T(0) ) - sust_(0,0,j,jp,0,0)  ) < 0.00000001 );

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
			{
				assert( std::abs( -( (j==jp)? j<2? T(1):T(-1) : T(0) )/4.0 - sust_(0,0,j,jp,0,0)  ) < 0.00000001 );
			}

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
			{
				assert( std::abs( -( (j==jp)? j%2==0? T(1):T(-1) : T(0) )/4.0 - sust_(0,0,j,jp,0,0)  ) < 0.00000001 );
			}

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
			{
				assert( std::abs( -( (j==jp)? j==1||j==2? T(-1):T(1) : T(0) )/4.0 - sust_(0,0,j,jp,0,0)  ) < 0.00000001 );
			}

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
			{
				assert( std::abs( -( (j==jp)? j==0? T(-1):T(1) : T(0) )/4.0 - sust_(0,0,j,jp,0,0)  ) < 0.00000001 );
			}

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
			{
				assert( std::abs( -( (j==jp)? j==1? T(-1):T(1) : T(0) )/4.0 - sust_(0,0,j,jp,0,0)  ) < 0.00000001 );
			}

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
								assert( std::abs( -( (j==jp)? v : T(0) )/4.0 - sust_(0,0,j,jp,l1*nOrb+l2,l3*nOrb+l4)  ) < 0.00000001 );
							}
	}

	const bT bandwidth = 20; // meV
	const bT mu = bandwidth*0.45;
	std::vector<size_t> spaceGrid = { 13, 11 };

	auto cos_bnd = [&] (std::vector<size_t> const& tuple,
			size_t l1, size_t a1, size_t s1,
			size_t l2, size_t a2, size_t s2){
		return ( a1 == a2 ? a1 == 0 ? 1.0 : -1.0 : 0.0 ) *  ( s1 == s2 ? 1.0 : 0.0 )
				*( 0.5*bandwidth*(std::cos( (2*M_PI*tuple[0])/spaceGrid[0] )
							+std::cos( (2*M_PI*tuple[1])/spaceGrid[1] ))-mu);
	};

	bT temperature = 1.0;
	GreensFunctionOrbital_test<T> gftest;
	gf = gftest.construct_gf_bnd(
			temperature, spaceGrid , 1,
			1,
			cos_bnd);

	gf.perform_space_fft();

	const bT kb = 0.86173324 ; // meV / K
	const bT beta = 1.0 / (kb * temperature);
	gf.transform_itime_Mfreq(beta);

	sust_.set_uninitialized();
	sust_.compute_from_gf( gf );

	sust_.perform_space_fft();

	std::ofstream file("/home/alinsch/codes/scallop/tests/sust_test.dat");
	for (size_t ik = 0 ; ik < sust_.get_spaceGrid_proc().get_num_k_grid(); ++ik)
	{
		auto tuple = sust_.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		file << tuple[0] << '\t' << tuple[1] << '\t' << sust_(ik,0,0,0,0,0)<< '\t' << sust_(ik,0,1,1,0,0)
				<< '\t' << sust_(ik,0,2,2,0,0)<< '\t' << sust_(ik,0,3,3,0,0) << '\n';
	}
	file.close();

	sust_.transform_itime_Mfreq(beta);
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
