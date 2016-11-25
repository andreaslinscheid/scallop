/*	This file SelfEnergy_test.cpp is part of scallop.
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
 *  Created on: Nov 24, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/SelfEnergy.h"
#include "scallop/gw_flex/GreensFunctionOrbital.h"
#include "scallop/gw_flex/GeneralizedSusceptibility.h"

namespace scallop
{
namespace gw_flex
{
namespace test
{

template<typename T>
class SelfEnergy_test
{
public:
	void test_all();

private:

	void creation_test();

};

template<typename T>
void SelfEnergy_test<T>::test_all()
{
	creation_test();
}

template<typename T>
void SelfEnergy_test<T>::creation_test()
{
	output::TerminalOut msg;
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	auto pauli_0 = [] (size_t i , size_t j) {return i == j ? T(1.0) 					: T(0.0);};
//	auto pauli_x = [] (size_t i , size_t j) {return i != j ? T(1.0) 					: T(0.0);};
	auto pauli_y = [] (size_t i , size_t j) {return i != j ? (i<j?T(0,1.0)	:T(0,-1.0))	: T(0.0);};
	auto pauli_z = [] (size_t i , size_t j) {return i == j ? (i==0?T(1.0)	:T(-1.0)) 	: T(0.0);};

	SelfEnergy<T> se;
	if ( mpi.get_nproc() == 1 )
	{
		std::vector<size_t> grid = {1,1};
		msg << "Testing basic construction of the self-energy at a single point in space and time";
		GreensFunctionOrbital<T> gf;
		typename auxillary::TemplateTypedefs<T>::scallop_vector data;
		gf.initialize(1,grid,1,true,false,data);

		msg << "\tcharge green's function:";
		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
						gf(0,0,0,a1,s1,0,a2,s2) = pauli_z(a1,a2) * pauli_0(s1,s2);

		GeneralizedSusceptibility<T> sust;
		sust.compute_from_gf( gf );

		se.set_to_zero();
		se.add_electronic_selfenergy( gf , sust);

		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
						assert( std::abs(( ( (a1==a2) && (s1==s2) )? a1 ==0? T(-1.0):T(1.0) : T(0) ) - se(0,0,0,a1,s1,0,a2,s2)  ) < 0.00000001 );

		msg << "\t singlet SC+charge+mz green's function (without recomputing the susceptibility):";
		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
						gf(0,0,0,a1,s1,0,a2,s2) = 1.0 * pauli_z(a1,a2) * pauli_0(s1,s2)
													-0.1*pauli_y(a1,a2) * pauli_y(s1,s2)
													+0.01*pauli_z(a1,a2) * pauli_z(s1,s2);

		se.set_to_zero();
		se.add_electronic_selfenergy( gf , sust);
		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
					{
						T val = ( (a1==a2) && (s1==s2) )? a1 ==0? T(-1.0):T(1.0) : T(0) ;
						val = ( (a1!=a2) && (s1!=s2) )? s2+a2 == s1+a1? T(0.05):T(-0.05) : val ;
						assert( std::abs(val - se(0,0,0,a1,s1,0,a2,s2)  ) < 0.00000001 );
					}

		size_t nO = 2;
		auto bnd = [] (size_t l1 , size_t l2 ){ return l1 == l2 ? l1 == 0 ? T(l1+1) : T(0,l1+1) : T(0);  };
		msg << "\t "<<nO<<" band charge green's function:";
		gf.initialize(1,grid,nO,true,false,data);
		for ( size_t l1 = 0 ; l1 < nO ; ++l1 )
			for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
				for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
					for ( size_t l2 = 0 ; l2 < nO ; ++l2 )
						for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
							for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
								gf(0,0,l1,a1,s1,l2,a2,s2) = bnd(l1,l2) * pauli_z(a1,a2) * pauli_0(s1,s2);

		//this sets the right dimensions
		sust.set_uninitialized();
		sust.compute_from_gf( gf );

		for ( size_t j = 0 ; j < 4 ; ++j )
			for ( size_t jp = 0 ; jp < 4 ; ++jp )
				for ( size_t l1 = 0 ; l1 < nO ; ++l1 )
					for ( size_t l2 = 0 ; l2 < nO ; ++l2 )
						for ( size_t l3 = 0 ; l3 < nO ; ++l3 )
							for ( size_t l4 = 0 ; l4 < nO ; ++l4 )
								sust(0,0,j,jp,l1,l2,l3,l4) = - 1.0 / 4.0 * ( l1 ==l3?1.0:0.0 ) * ( l2 ==l4?1.0:0.0 ) * ( j ==jp?1.0:0.0 );

		se.set_uninitialized();
		se.add_electronic_selfenergy( gf , sust);

		for ( size_t l1 = 0 ; l1 < nO ; ++l1 )
			for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
				for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
					for ( size_t l2 = 0 ; l2 < nO ; ++l2 )
						for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
							for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
							{
								T val = ( l1 == l2?1.0:0.0 );
								val *= -1.0*(gf(0,0,0,a1,s1,0,a2,s2)+gf(0,0,1,a1,s1,1,a2,s2));
								assert( std::abs(val - se(0,0,l1,a1,s1,l2,a2,s2)  ) < 0.00000001 );
							}
	}

}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
