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
#include "scallop/output/TerminalOut.h"
#include "scallop/parallel/MPIModule.h"
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

	void one_loop_z_test();

	void interpolate_test();

};

template<typename T>
void SelfEnergy_test<T>::test_all()
{
	creation_test();
	interpolate_test();
	one_loop_z_test();
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

		SpinSusceptibility<T> sust;
		sust.compute_from_gf( gf );

		se.add_electronic_selfenergy( gf , sust);

		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
						assert( std::abs(( ( (a1==a2) && (s1==s2) )? a1 ==0? T(-1.0):T(1.0) : T(0) ) + se(0,0,0,a1,s1,0,a2,s2)  ) < 0.00000001 );

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
						val = -1.0*(( (a1!=a2) && (s1!=s2) )? s2+a2 == s1+a1? T(0.05):T(-0.05) : val );
						assert( std::abs(val - se(0,0,0,a1,s1,0,a2,s2)  ) < 0.00000001 );
					}

		msg << "\tFrequency part of the green's function:";
		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
						gf(0,0,0,a1,s1,0,a2,s2) = pauli_0(a1,a2) * pauli_0(s1,s2);
		se.set_to_zero();
		se.add_electronic_selfenergy( gf , sust);
		for ( size_t a1 = 0 ; a1 < 2 ; ++a1 )
			for ( size_t s1 = 0 ; s1 < 2 ; ++s1 )
				for ( size_t a2 = 0 ; a2 < 2 ; ++a2 )
					for ( size_t s2 = 0 ; s2 < 2 ; ++s2 )
					{
						T val = -1.0*(( (a1==a2) && (s1==s2) )? -T(1.0) : T(0) );
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
								val *= (gf(0,0,0,a1,s1,0,a2,s2)+gf(0,0,1,a1,s1,1,a2,s2));
								assert( std::abs(val - se(0,0,l1,a1,s1,l2,a2,s2)  ) < 0.00000001 );
							}
	}

}

template<typename T>
void SelfEnergy_test<T>::one_loop_z_test()
{
	output::TerminalOut msg;
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	msg << "Testing the frequency renormalzation part of a one band cosine model.";

	auxillary::BasicFunctions bf;
	auto scallopPath = bf.get_scallop_path();
	std::string const filenameModel = scallopPath+"examples/oneBndCosine/model_test.dat";
	std::string const filenameInter = scallopPath+"examples/oneBndCosine/input_data/Isf.dat";
	if ( mpi.ioproc() )
	{
		std::ofstream file( filenameModel.c_str() );
		file << "model OneBandCosine\nt=290\ntp=-37.5";
		file.close();
		file.open( filenameInter.c_str() );
		file << "1 4\n"
				"0\t0\t0\t0\t0\t0\t200.0\t0.0\n"
				"1\t1\t0\t0\t0\t0\t200.0\t0.0\n"
				"2\t2\t0\t0\t0\t0\t200.0\t0.0\n"
				"3\t3\t0\t0\t0\t0\t200.0\t0.0\n";
		file.close();
	}

	typedef typename auxillary::TypeMapComplex<T>::type bT;
	const bT temperature = 10;
	const bT beta = 1.0 / (auxillary::Constants<bT>::kBoltzmannMEV * temperature);
	const size_t nM = 512;

	KohnShamBandStructure<T> ksBnd;
	std::vector<size_t> grid = {32, 32};
	ksBnd.initialize_from_file( grid, filenameModel );
	ksBnd.adjust_filling(0.875,beta);

	msg << "\tFilling is " << ksBnd.compute_N_electrons( beta ) << " and temperature " << temperature << "K";

	KohnShamGreensFunctionOrbital<T> ksGF;
	ksGF.set_from_KS_bandstructure(true,nM,beta,ksBnd);
	ksGF.perform_space_fft( );

	InteractionMatrix<T> Isf;
	Isf.init_file( filenameInter );

	SpinSusceptibility<T> suscSF,enhSF;
	suscSF.compute_from_gf( ksGF );
	suscSF.transform_itime_Mfreq( beta );
	suscSF.perform_space_fft( );
	enhSF = suscSF;
	suscSF.spin_RPA_enhancement( Isf );
	enhSF.spin_RPA_enhancement( Isf, true );
	std::ofstream fileBla(scallopPath+"examples/oneBndCosine/output_data/bla.dat");
	for (size_t ik = 0 ; ik < suscSF.get_spaceGrid_proc().get_num_k_grid(); ++ik)
	{
		auto tuple = suscSF.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		fileBla << bT(tuple[0])/bT(grid[0]) << '\t' << bT(tuple[1])/bT(grid[0]) << '\t' << std::real(enhSF(ik,0,0,0,0,0))
			<< '\t' << std::real(enhSF(ik,0,1,1,0,0))
			<< '\t' << std::real(enhSF(ik,0,2,2,0,0))
			<< '\t' << std::real(enhSF(ik,0,3,3,0,0))<< std::endl;
		if ( tuple[1]==grid[1]-1 )
			fileBla  << std::endl;
	}
	fileBla.close();
	suscSF.transform_itime_Mfreq( beta );
	suscSF.perform_space_fft( );

	SelfEnergy<T> se;
	se.add_electronic_selfenergy(ksGF,suscSF);
	se.transform_itime_Mfreq( beta );
	se.perform_space_fft( );

	assert( se.is_in_k_space() );
	assert( not se.is_in_time_space() );

	//Check Gamma and pi 0
	size_t mid = grid[1]*(grid[0]/2)+0;
	std::ofstream fileGamma(scallopPath+"examples/oneBndCosine/output_data/z_gamma.dat");
	std::ofstream filePiPi(scallopPath+"examples/oneBndCosine/output_data/z_pi.dat");
	for (size_t ik = 0 ; ik < se.get_spaceGrid_proc().get_num_k_grid(); ++ik)
	{
		auto tuple = se.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ictotal = se.get_spaceGrid_proc().k_xyz_to_conseq( tuple );
		if ( ictotal == 0 )
			for (size_t iw = 0 ; iw < nM; ++iw)
			{
				auto ptr = se.read_phs_grid_ptr_block(ik,iw);
				T z = T(0);
				for (size_t ns = 0 ; ns < 4; ++ns)
					z += 0.25*ptr[ns*4+ns];
				z = bT(1.0)+T(0,1.0)*std::imag(z)/auxillary::BasicFunctions::matzubara_frequency_of_index(iw,nM,beta);
				fileGamma << std::imag(auxillary::BasicFunctions::matzubara_frequency_of_index(iw,nM,beta)) << '\t'
						<< std::real(z) <<std::endl;
			}
		if ( ictotal == mid )
			for (size_t iw = 0 ; iw < nM; ++iw)
			{
				auto ptr = se.read_phs_grid_ptr_block(ik,iw);
				T z = T(0);
				for (size_t ns = 0 ; ns < 4; ++ns)
					z += 0.25*ptr[ns*4+ns];
				z = bT(1.0)+T(0,1.0)*std::imag(z)/auxillary::BasicFunctions::matzubara_frequency_of_index(iw,nM,beta);
				filePiPi << std::imag(auxillary::BasicFunctions::matzubara_frequency_of_index(iw,nM,beta)) << '\t'
						<< std::real(z) <<std::endl;
			}
	}
	fileGamma.close();
	filePiPi.close();

}

template<typename T>
void SelfEnergy_test<T>::interpolate_test()
{
	output::TerminalOut msg;
	msg << "Testing frequency interpolation of the self-energy";

	parallel::MPIModule const & mpi = parallel::MPIModule::get_instance();

	if ( mpi.get_nproc() == 1 )
	{
		size_t const nM = 10;
		auto beta = auxillary::BasicFunctions::inverse_temperature( 10.0 );
		auto firstGrid = auxillary::BasicFunctions::matzubara_frequency_array(nM,beta);
		std::vector<size_t> grid = { 1, 1 };
		typename auxillary::TemplateTypedefs<T>::scallop_vector data(nM*4*4);
		for ( size_t n = 0 ; n < nM ; ++n)
		{
			int i = n < nM/2 ? n : static_cast<int>(n)-static_cast<int>(nM);
			data[n*16] = i;
		}
		SelfEnergy<T> se;
		se.initialize(
				nM,
				grid,
				1,
				false,
				true,
				data);

		auto beta2 = auxillary::BasicFunctions::inverse_temperature( 5.0 );
		size_t const nM2 = 20;
		auto secondGrid = auxillary::BasicFunctions::matzubara_frequency_array(nM2,beta2);
		se.linear_interpolate_frequency(firstGrid,secondGrid);

		//Check if this is a straight line, at least in the inner range
		typedef typename auxillary::TypeMapComplex<T>::type bT;
		for ( size_t n = 2 ; n < nM2-2 ;++n)
		{
			size_t r = n < nM2/2 ?
					static_cast<int>(n)+nM2/2
					:static_cast<int>(n)-static_cast<int>(nM2/2);
			bT ref = nM/2*(bT(n)-bT(nM2/2))/bT(nM2/2)-0.25;
			if ( std::abs( ref-std::real(se(0,r,0,0)) ) > 0.00001 )
			{
				//error_handling::Error("Check failed, linear interpolation self-energy");
				std::cout << ref << '\t' << std::real(se(0,r,0,0)) <<std::endl;
			}
		}
	}
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
