/*	This file FFTBase_test.hpp is part of scallop.
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
 *  Created on: Nov 10, 2016
 *      Author: A. Linscheid
 */

#include "FFTBase_test.h"
#include "scallop/auxillary/TemplateTypedefs.h"
#include <fstream>

namespace scallop
{
namespace gw_flex
{
namespace test
{

template<typename T>
void FFTBase_test<T>::test_all()
{
	test_initialize();
	test_copy_assign_move();
	test_fft_space();
}

template<typename T>
void FFTBase_test<T>::test_copy_assign_move()
{
	output::TerminalOut msg;
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	msg << "Testing copy, move and assign";

	std::vector<size_t> spaceGridTotal = { 4, 3 };
	parallel::GridDistribution<T> gd;
	gd.distribute_grid(spaceGridTotal);
	bool dataIsInKSpace = true;
	bool dataIsInTimeSpace= true;
	size_t dimTimeFT = 1;
	size_t blockSize =1 ;
	typename auxillary::TemplateTypedefs<T>::scallop_vector data( gd.get_num_k_grid()*dimTimeFT*blockSize ,T(1));
	int i = 0;
	for ( auto &d : data )
		d = T(++i,mpi.get_mpi_me());

	FFTBase<T> ffttest;
	ffttest.initialize( data, dataIsInKSpace, dataIsInTimeSpace, spaceGridTotal, dimTimeFT, blockSize);

	FFTBase<T> ffttestCpy(ffttest);
	T diff = 0;
	for ( size_t i = 0 ; i < data.size(); ++i)
		diff += abs(ffttest.read_data_ptr_block(0,0)[i]-ffttestCpy.read_data_ptr_block(0,0)[i]);

	assert( (diff.real() == 0) && (diff.imag() == 0) );

	FFTBase<T> ffttestAsgn = ffttest;
	diff = 0;
	for ( size_t i = 0 ; i < data.size(); ++i)
		diff += abs(ffttest.read_data_ptr_block(0,0)[i]-ffttestAsgn.read_data_ptr_block(0,0)[i]);

	assert( (diff.real() == 0) && (diff.imag() == 0) );

	FFTBase<T> ffttestMove = std::move(ffttest);
	diff = 0;
	for ( size_t i = 0 ; i < data.size(); ++i)
		diff += abs(ffttestMove.read_data_ptr_block(0,0)[i]-ffttestAsgn.read_data_ptr_block(0,0)[i]);

	assert( (diff.real() == 0) && (diff.imag() == 0) );
}

template<typename T>
void FFTBase_test<T>::test_initialize()
{
	output::TerminalOut msg;
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	if ( mpi.get_nproc() < 2 )
	{
		msg << "Testing FFTBase initialization for 1 grid point";
		FFTBase<T> ffttest();

		typename auxillary::TemplateTypedefs<T>::scallop_vector data(1,T(1));
		bool dataIsInKSpace = true;
		bool dataIsInTimeSpace= true;
		std::vector<size_t> spaceGridTotal = { 1, 1 };
		size_t dimTimeFT = 1;
		size_t blockSize =1 ;

		FFTBase<T> ffttest2;
		ffttest2.initialize( data, dataIsInKSpace, dataIsInTimeSpace, spaceGridTotal, dimTimeFT, blockSize);

		typename auxillary::TemplateTypedefs<T>::scallop_vector data2(16,T(1));
		blockSize = 0;
		for ( auto &d : data2 )
			d = ++blockSize;
		ffttest2.initialize( data2, dataIsInKSpace, dataIsInTimeSpace, spaceGridTotal, dimTimeFT, blockSize);
		T diff = 0;
		for ( size_t i = 0 ; i < blockSize; ++i)
			diff += T(std::real(ffttest2.read_data_ptr_block(0,0)[i]-T(i+1)),
					std::imag(ffttest2.read_data_ptr_block(0,0)[i]-T(i+1)));

		assert( (diff.real() == 0) && (diff.imag() == 0) );
	}

	msg << "Testing FFTBase initialization a 4x3 grid";

	std::vector<size_t> spaceGridTotal = { 4, 3 };
	parallel::GridDistribution<T> gd;
	gd.distribute_grid(spaceGridTotal);
	bool dataIsInKSpace = true;
	bool dataIsInTimeSpace= true;
	size_t dimTimeFT = 1;
	size_t blockSize =1 ;
	typename auxillary::TemplateTypedefs<T>::scallop_vector data( gd.get_num_k_grid()*dimTimeFT*blockSize ,T(1));
	int i = 0;
	for ( auto &d : data )
		d = T(++i,mpi.get_mpi_me());

	FFTBase<T> ffttest;
	ffttest.initialize( data, dataIsInKSpace, dataIsInTimeSpace, spaceGridTotal, dimTimeFT, blockSize);

	T diff = 0;
	for ( size_t i = 0 ; i < data.size(); ++i)
	{
		diff += T(std::real(ffttest.read_phs_grid_ptr_block(i,0)[0]-T(i+1,mpi.get_mpi_me())),
				std::imag(ffttest.read_phs_grid_ptr_block(i,0)[0]-T(i+1,mpi.get_mpi_me())));
	}

	mpi.sum( diff );
	msg << "Difference in the data passed to the initialization routine and the set values: " << diff;

	mpi.barrier();
	assert( (diff.real() == 0) && (diff.imag() == 0) );

	FFTBase<T> ffttestEmpty;
	typename auxillary::TemplateTypedefs<T>::scallop_vector dataEmpty;
	ffttest.initialize( dataEmpty, dataIsInKSpace, dataIsInTimeSpace, spaceGridTotal, dimTimeFT, blockSize);
	diff = 0;
	for ( size_t i = 0 ; i < data.size(); ++i)
		diff += T(std::real(ffttest.read_phs_grid_ptr_block(i,0)[0]-T(0)),
				std::imag(ffttest.read_phs_grid_ptr_block(i,0)[0]-T(0)));

	assert( (diff.real() == 0) && (diff.imag() == 0) );
}

template<typename T>
void FFTBase_test<T>::test_fft_space()
{
	output::TerminalOut msg;
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	msg << "Testing FFTBase space Fourier transform. Currently using "
			<< mpi.get_nproc() << " processors.";
	mpi.barrier();

	std::vector<size_t> testgrid = {10, 11};
	parallel::GridDistribution<T> gd;
	gd.distribute_grid(testgrid);

	FFTBase<T> ffttest;
	msg << "Testing transform of a cosine kx plus cosine ky in k space to lattice space";
	const size_t nM = 1;
	const bT bandwidth = 20; // meV
	const bT mu = bandwidth*0.9;
	auto cos_test = [&] (size_t ikx, size_t iky)
		{
			return 0.5*bandwidth*(std::cos( (2*M_PI*ikx)/testgrid[0] )
				+std::cos( (2*M_PI*iky)/testgrid[1] ))-mu;
		};

	typename auxillary::TemplateTypedefs<T>::scallop_vector data( gd.get_num_k_grid() );
	for (size_t ik = 0 ; ik < gd.get_num_k_grid() ; ++ik)
	{
		auto tuple= gd.k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple.front();
		size_t iky = tuple.back();
		data[ik] = cos_test(ikx,iky);
	}

	ffttest.initialize( data, true, true, testgrid, nM, 1 );
	size_t dRx = ffttest.get_spaceGrid_proc().get_grid()[0];
	size_t dRy = ffttest.get_spaceGrid_proc().get_grid()[1];

	for ( size_t iproc = 0 ; iproc < mpi.get_nproc(); ++iproc )
	{
		if (mpi.get_mpi_me() == iproc)
		{
			std::cout << "Processor "<< mpi.get_mpi_me() << " k grid dimension: \n";
			for (auto k : ffttest.get_spaceGrid_proc().get_k_grid() )
				std::cout << k << "\t";
			std::cout << '\n' ;
		}
		mpi.barrier();
	}

	ffttest.perform_space_fft();

	T diff = 0;
	for (size_t iR = 0 ; iR < ffttest.get_spaceGrid_proc().get_num_R_grid(); ++iR)
	{
		auto tuple= ffttest.get_spaceGrid_proc().R_conseq_local_to_xyz_total( iR );
		size_t iRx = tuple.front();
		size_t iRy = tuple.back();

		T analytic = T(0);
		if ((iRx==iRy)&&(iRx==0))
			analytic = -mu;

		if ( (iRx == 0)&&((iRy == 1)||(iRy == dRy-1)) )
			analytic = 0.25*bandwidth;
		if ( (iRy == 0)&&((iRx == 1)||(iRx == dRx-1)) )
			analytic = 0.25*bandwidth;

		diff += std::abs(std::real(ffttest.read_phs_grid_ptr_block(iR,0)[0])-std::real(analytic));
		diff += std::abs(std::imag(ffttest.read_phs_grid_ptr_block(iR,0)[0])-std::imag(analytic));
	}
	mpi.sum( diff );
	msg << "Difference between the analytic results and the numerical routine: " << diff;

	mpi.barrier();
	assert( (diff.real() < 0.0000001) && (diff.imag() < 0.0000001) );

	msg << "Transforming back to k space and comparing ...";
	ffttest.perform_space_fft();
	diff = 0;

	for (size_t ik = 0 ; ik < ffttest.get_spaceGrid_proc().get_num_k_grid(); ++ik)
	{
		auto tuple= ffttest.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple.front();
		size_t iky = tuple.back();

		T analytic = cos_test(ikx,iky);

		diff += std::abs(std::real(ffttest.read_phs_grid_ptr_block(ik,0)[0])-std::real(analytic));
		diff += std::abs(std::imag(ffttest.read_phs_grid_ptr_block(ik,0)[0])-std::imag(analytic));
	}

	mpi.sum( diff );
	msg << "Difference between the analytic results and the numerical routine: " << diff;

	mpi.barrier();
	assert( (diff.real() < 0.0000001) && (diff.imag() < 0.0000001) );

	msg << "Testing transform of a cosine kx plus cosine ky band each in one of two components in k space to lattice space";
	data = typename auxillary::TemplateTypedefs<T>::scallop_vector( gd.get_num_k_grid()*2 );
	for (size_t ik = 0 ; ik < gd.get_num_k_grid() ; ++ik)
	{
		auto tuple= gd.k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple.front();
		size_t iky = tuple.back();
		data[ik*2 + 0] = cos_test(ikx,iky);
		data[ik*2 + 1] = cos_test(ikx,iky);
	}
	ffttest.initialize( data, true, true, testgrid, nM, 2 );

	ffttest.perform_space_fft();

	diff = 0;
	for (size_t iR = 0 ; iR < ffttest.get_spaceGrid_proc().get_num_R_grid(); ++iR)
	{
		auto tuple= ffttest.get_spaceGrid_proc().R_conseq_local_to_xyz_total( iR );
		size_t iRx = tuple.front();
		size_t iRy = tuple.back();

		T analytic = T(0);
		if ((iRx==iRy)&&(iRx==0))
			analytic = -mu;

		if ( (iRx == 0)&&((iRy == 1)||(iRy == dRy-1)) )
			analytic = 0.25*bandwidth;
		if ( (iRy == 0)&&((iRx == 1)||(iRx == dRx-1)) )
			analytic = 0.25*bandwidth;

		diff += std::abs(std::real(ffttest.read_phs_grid_ptr_block(iR,0)[0])-std::real(analytic));
		diff += std::abs(std::imag(ffttest.read_phs_grid_ptr_block(iR,0)[0])-std::imag(analytic));

		diff += std::abs(std::real(ffttest.read_phs_grid_ptr_block(iR,0)[1])-std::real(analytic));
		diff += std::abs(std::imag(ffttest.read_phs_grid_ptr_block(iR,0)[1])-std::imag(analytic));
	}

	mpi.sum( diff );
	msg << "Difference between the analytic results and the numerical routine: " << diff;
	mpi.barrier();
	assert( (diff.real() < 0.0000001) && (diff.imag() < 0.0000001) );


	msg << "Doing the same thing as before, but now the two components are in the frequency dimension";
	ffttest.initialize( data, true, true, testgrid, 2, 1 );

	ffttest.perform_space_fft();

	diff = 0;
	for (size_t iR = 0 ; iR < ffttest.get_spaceGrid_proc().get_num_R_grid(); ++iR)
	{
		auto tuple= ffttest.get_spaceGrid_proc().R_conseq_local_to_xyz_total( iR );
		size_t iRx = tuple.front();
		size_t iRy = tuple.back();

		T analytic = T(0);
		if ((iRx==iRy)&&(iRx==0))
			analytic = -mu;

		if ( (iRx == 0)&&((iRy == 1)||(iRy == dRy-1)) )
			analytic = 0.25*bandwidth;
		if ( (iRy == 0)&&((iRx == 1)||(iRx == dRx-1)) )
			analytic = 0.25*bandwidth;

		diff += std::abs(std::real(ffttest.read_phs_grid_ptr_block(iR,0)[0])-std::real(analytic));
		diff += std::abs(std::imag(ffttest.read_phs_grid_ptr_block(iR,0)[0])-std::imag(analytic));

		diff += std::abs(std::real(ffttest.read_phs_grid_ptr_block(iR,1)[0])-std::real(analytic));
		diff += std::abs(std::imag(ffttest.read_phs_grid_ptr_block(iR,1)[0])-std::imag(analytic));

	}

	msg << "Difference between the analytic results and the numerical routine: " << diff;
	assert( (diff.real() < 0.0000001) && (diff.imag() < 0.0000001) );

	msg << "Transforming back ...";
	ffttest.perform_space_fft();
	diff = 0;
	for (size_t ik = 0 ; ik < ffttest.get_spaceGrid_proc().get_num_k_grid(); ++ik)
	{
		auto tuple= ffttest.get_spaceGrid_proc().k_conseq_local_to_xyz_total( ik );
		size_t ikx = tuple.front();
		size_t iky = tuple.back();

		T analytic = 0.5*bandwidth*(std::cos( (2*M_PI*ikx)/dRx )
			+std::cos( (2*M_PI*iky)/dRy ))-mu;

		diff += std::abs(std::real(ffttest.read_phs_grid_ptr_block(ik,0)[0])-std::real(analytic));
		diff += std::abs(std::imag(ffttest.read_phs_grid_ptr_block(ik,0)[0])-std::imag(analytic));

		diff += std::abs(std::real(ffttest.read_phs_grid_ptr_block(ik,1)[0])-std::real(analytic));
		diff += std::abs(std::imag(ffttest.read_phs_grid_ptr_block(ik,1)[0])-std::imag(analytic));
	}

	msg << "Difference between the analytic results and the numerical routine: " << diff;
	assert( (diff.real() < 0.0000001) && (diff.imag() < 0.0000001) );

}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
