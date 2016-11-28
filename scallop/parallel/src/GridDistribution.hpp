/*	This file GridDistribution.hpp is part of scallop.
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
 *  Created on: Nov 4, 2016
 *      Author: A. Linscheid
 */

#include "scallop/parallel/GridDistribution.h"
#include "scallop/error_handling/Error.h"
#include "scallop/parallel/MPIModule.h"
#include "scallop/output/TerminalOut.h"
#include <set>

namespace scallop
{
namespace parallel
{

template<typename T>
GridDistribution<T>::GridDistribution()
{

}

template<typename T>
void GridDistribution<T>::distribute_grid( std::vector<size_t> totalGrid )
{
	totalGrid_ = std::move(totalGrid);
	size_t lastD = totalGrid_.size()-1;

	if ( totalGrid_.size() < 2 )
		error_handling::Error( "Cannot distribute grids smaller than 2 D" );

	//prepare the k grid per processor
	this->distribute_dim( totalGrid_[0], parallelMapk_ );

	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	size_t dimKxProc = static_cast<size_t>(
			static_cast<int>( parallelMapk_[ mpi.get_mpi_me() ].second)
			- static_cast<int>( parallelMapk_[ mpi.get_mpi_me() ].first ) );
	procGridk_ = totalGrid_;
	procGridk_[0] = dimKxProc;

	//prepare the R grid per processor
	this->distribute_dim( totalGrid_[lastD], parallelMapR_ );

	size_t dimRLProc = static_cast<size_t>(
			static_cast<int>( parallelMapR_[ mpi.get_mpi_me() ].second)
			- static_cast<int>( parallelMapR_[ mpi.get_mpi_me() ].first ) );
	procGridR_ = totalGrid_;
	procGridR_[lastD] = dimRLProc;

	nK_ = 1;
	for ( auto k : procGridk_ )
		nK_ *= k;

	nR_ = 1;
	for ( auto R : procGridR_ )
		nR_ *= R;

	nGridTotal_ = 1;
	for ( auto kR : totalGrid_ )
		nGridTotal_ *= kR;

	//Determine the data grid. This is the minimal grid which
	// 1. has the same number of elements both in R and k and
	// 2. has the same size on all processors
	size_t nP = mpi.get_nproc();
	size_t nfT = totalGrid_[0];
	size_t nlT = totalGrid_[lastD];
	procGridDatak_ = totalGrid_;
	procGridDataR_ = totalGrid_;
	this->tuple_transpose(procGridDataR_);
	procGridDatak_[0] = ( nfT % nP == 0 ) ? nfT / nP : nfT / nP + 1 ;
	procGridDatak_[lastD] = ( nlT % nP == 0 ) ? nlT  : nlT + 1 ;
	procGridDataR_[0] = ( nfT % nP == 0 ) ? nfT : nfT + 1 ;
	procGridDataR_[lastD] = ( nlT % nP == 0 ) ? nlT / nP : nlT / nP + 1 ;
	nProcGridData_ = 1;
	for ( auto kRD : procGridDatak_ )
		nProcGridData_ *= kRD;

	size_t test = 1;
	for ( auto kRD : procGridDataR_ )
		test *= kRD;
	if ( nProcGridData_ != test )
		error_handling::Error("It seems R and k grid do not have the same number of data elements. The algorithm is bad.");
}

template<typename T>
void GridDistribution<T>::distribute_dim(
		size_t N,
		std::vector< std::pair<size_t,size_t> > & procDimDistribution ) const
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	size_t nP = mpi.get_nproc();

	size_t dimMe = N / nP;

	//assign the N left over grid points to the first N processors
	// so that each one will handle 1 more than the default dimMe at most
	size_t leftover = static_cast<size_t>(
			static_cast<int>(N) - static_cast<int>( dimMe*nP ) );

	procDimDistribution = std::vector<std::pair<size_t,size_t> >( nP );
	for ( size_t i = 0 ; i < nP ; i++  )
	{
		size_t addtoDefaultDim = i < leftover ? i : leftover;
		size_t procIGridStart = N / nP * i + addtoDefaultDim;
		procDimDistribution[i].first = procIGridStart;
		procDimDistribution[i].second = procIGridStart + dimMe + (i < leftover ? 1 : 0);
	}
}

template<typename T>
size_t GridDistribution<T>::get_dim() const
{
	return totalGrid_.size();
}

template<typename T>
std::vector<size_t> const& GridDistribution<T>::get_grid() const
{
	return totalGrid_;
}

template<typename T>
size_t  GridDistribution<T>::get_num_grid() const
{
	return nGridTotal_;
}

template<typename T>
size_t GridDistribution<T>::get_num_grid_data() const
{
	return nProcGridData_;
}

template<typename T>
std::vector<size_t> const& GridDistribution<T>::get_k_grid() const
{
	return procGridk_;
}

template<typename T>
size_t GridDistribution<T>::get_num_k_grid() const
{
	return nK_;
}

template<typename T>
std::vector<size_t> const& GridDistribution<T>::get_R_grid() const
{
	return procGridR_;
}

template<typename T>
size_t GridDistribution<T>::get_num_R_grid() const
{
	return nR_;
}

template<typename T>
size_t GridDistribution<T>::k_conseq_local_to_data_conseq( size_t ik ) const
{
	this->general_conseq_to_xyz_column_major( procGridk_, tupleBuff_, ik );
	this->general_xyz_to_conseq_column_major( procGridDatak_, tupleBuff_, ik );
	return ik;
}

template<typename T>
size_t GridDistribution<T>::R_conseq_local_to_data_conseq( size_t iR ) const
{
	this->general_conseq_to_xyz_row_major( procGridR_, tupleBuff_, iR );
	this->general_xyz_to_conseq_row_major( procGridDataR_, tupleBuff_, iR );
	return iR;
}

template<typename T>
std::vector<size_t> & GridDistribution<T>::k_conseq_to_xyz( size_t ik ) const
{
	this->general_conseq_to_xyz_column_major( procGridk_, tupleBuff_, ik );
	return tupleBuff_;
}

template<typename T>
std::vector<size_t> & GridDistribution<T>::R_conseq_to_xyz( size_t iR ) const
{
	this->general_conseq_to_xyz_row_major( procGridR_, tupleBuff_, iR);
	return tupleBuff_;
}

template<typename T>
size_t GridDistribution<T>::R_xyz_total_to_conseq_local( std::vector<size_t> tuple ) const
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	int coord = static_cast<int>(tuple.back())-static_cast<int>(parallelMapR_[ mpi.get_mpi_me() ].first);
#ifdef DEBUG_BUILD
	if ( (coord < 0) || (coord >= static_cast<int>(parallelMapR_[ mpi.get_mpi_me() ].second)) )
		error_handling::Error("The tuple passed corresponds to a grid vector out of range of this processor.");
#endif
	tuple.back() = static_cast<size_t>(coord);
	return this->R_xyz_to_conseq( tuple );
}

template<typename T>
size_t GridDistribution<T>::R_xyz_to_conseq( std::vector<size_t> const& tuple ) const
{
	size_t index;
	this->general_xyz_to_conseq_row_major( procGridR_, tuple, index );
	return index;
}

template<typename T>
size_t GridDistribution<T>::k_xyz_to_conseq( std::vector<size_t> const& tuple ) const
{
	size_t index;
	this->general_xyz_to_conseq_column_major( procGridk_, tuple, index );
	return index;
}

template<typename T>
std::vector<size_t> & GridDistribution<T>::k_conseq_local_to_xyz_total( size_t ik ) const
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	tupleBuff_ = this->k_conseq_to_xyz( ik );
	tupleBuff_.front() += parallelMapk_[ mpi.get_mpi_me() ].first;
	return tupleBuff_;
}

template<typename T>
std::vector<size_t> & GridDistribution<T>::R_conseq_local_to_xyz_total( size_t iR ) const
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	tupleBuff_ = this->R_conseq_to_xyz( iR );
	tupleBuff_.back() += parallelMapR_[ mpi.get_mpi_me() ].first;
	return tupleBuff_;
}

template<typename T>
size_t GridDistribution<T>::get_num_grid_total() const
{
	return nGridTotal_;
}

template<typename T>
size_t GridDistribution<T>::get_inverse_index_k( size_t ik ) const
{
	tupleBuff_ = this->k_conseq_local_to_xyz_total(ik);
	for ( size_t i = 0; i < tupleBuff_.size(); ++i )
		tupleBuff_[i] = static_cast<size_t>(static_cast<int>(totalGrid_[i])-static_cast<int>(tupleBuff_[i]));
	return this->k_xyz_to_conseq( tupleBuff_ );
}

template<typename T>
size_t GridDistribution<T>::get_inverse_index_R( size_t iR ) const
{
	tupleBuff_ = this->R_conseq_local_to_xyz_total(iR);
	for ( size_t i = 0; i < tupleBuff_.size(); ++i )
		tupleBuff_[i] = static_cast<size_t>(static_cast<int>(totalGrid_[i])-static_cast<int>(tupleBuff_[i]));
	return this->R_xyz_to_conseq( tupleBuff_ );
}

template<typename T>
void GridDistribution<T>::grid_data_transposition(
		bool startFromKSpace,
		typename auxillary::TemplateTypedefs<T>::scallop_vector & data,
		size_t blockSize) const
{
	//First we need to rearrange the ordering so that the parallelized grid
	//dimension is the fastest running (the corresponding data blocks are adjacent) on each processor

	//sort to the new alignment
	//The following map tells for every entry the index that is supposed to overwrite the data at the present one.
	//Thus, we can buffer data at one index, and then walk though the array and always overwrite the present data
	//	with the next entry. In the end, because the circle must close, we can overwrite the last index with the
	//	buffer.
	std::vector<size_t> oldToNewAligned(nProcGridData_);
	std::vector<size_t> newToOldAligned(nProcGridData_);
	for ( size_t ig = 0; ig < nProcGridData_; ++ig)
	{
		size_t newAlignedIndex = 0;
		if ( startFromKSpace )
		{
			this->general_conseq_to_xyz_row_major(procGridDatak_,tupleBuff_,ig);
			this->general_xyz_to_conseq_column_major(procGridDatak_,tupleBuff_,newAlignedIndex);
		}
		else
		{
			this->general_conseq_to_xyz_column_major(procGridDataR_,tupleBuff_,ig);
			this->general_xyz_to_conseq_row_major(procGridDataR_,tupleBuff_,newAlignedIndex);
		}
		oldToNewAligned[ig] = newAlignedIndex;
	}

	this->redistribute_locally( data, oldToNewAligned, blockSize );

#ifdef MPI_PARALLEL
	//Now we exchange the data among processors
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	size_t procBlockSize = nProcGridData_*blockSize/mpi.get_nproc();
	mpi.all_to_all<T>( data, procBlockSize );

	//At this point, the data is incoherent. At each processor we have NProc patches of sorted data
	//The ordering of data within such a patch is described by blockLocalGrid, the local grid
	//which is a minimum of all data grids.
	std::vector<size_t> blockLocalGrid;
	size_t orderedDim;
	if ( startFromKSpace )
	{
		blockLocalGrid = procGridDatak_;
		blockLocalGrid.back() = procGridDataR_.back();
		orderedDim = blockLocalGrid.front();
	}
	else
	{
		blockLocalGrid = procGridDataR_;
		blockLocalGrid.front() = procGridDatak_.front();
		orderedDim = blockLocalGrid.back();
	}

	size_t nPerBlock = 1;
	for ( auto i : blockLocalGrid)
		nPerBlock *= i;
	nPerBlock /= orderedDim;

	std::vector<size_t> mapSyncToNewLayout(nProcGridData_);
	size_t count=0;
	for ( size_t ib = 0 ; ib < nPerBlock; ++ib )
		for ( size_t ip = 0; ip < mpi.get_nproc(); ++ip)
			for ( size_t ig = 0 ; ig < orderedDim; ++ig )
			{
				mapSyncToNewLayout[count++]=(ip*nPerBlock+ib)*orderedDim+ig;
			}

	this->redistribute_locally( data, mapSyncToNewLayout, blockSize );

#endif
}

template<typename T>
void GridDistribution<T>::tuple_transpose( std::vector<size_t> & tuple ) const
{
	int dim = static_cast<int>(tuple.size());
	std::vector<size_t> tmp = tuple;
	for ( int i = 0; i < dim ; ++i )
		tuple[i] = tmp[dim-1-i];
}

template<typename T>
void GridDistribution<T>::general_xyz_to_conseq_row_major(
		std::vector<size_t> const& grid,
		std::vector<size_t> const& indexTouple,
		size_t & conseq_index) const
{
#ifdef DEBUG_BUILD
	if (indexTouple.size() != grid.size())
		error_handling::Error("Incorrect in put, dimension must match!");
#endif
	conseq_index = indexTouple.back();
	for ( int i = static_cast<int>(indexTouple.size())-2 ; i >= 0; i--)
	{
		conseq_index *= grid[i];
		conseq_index += indexTouple[i];
	}
}

template<typename T>
void GridDistribution<T>::general_conseq_to_xyz_row_major(
		std::vector<size_t> const& grid,
		std::vector<size_t> & indexTouple,
		size_t conseq_index) const
{
	//reverse-engineer expressions like
	//	index = x + dx*( y + dy*( z ) )
	//	      = x + dx*y + dx*dy*z
	//for x,y,z using integer division

	int d = static_cast<int>(grid.size());
	std::vector<int> dimProd( d , 1 );
	for (int i = 1 ; i < d ; ++i)
		dimProd[i] = grid[ i-1 ]*dimProd[i-1];

	indexTouple = std::vector<size_t>( d , conseq_index );
	for ( int i = d-1 ; i >= 1; --i )
	{
		indexTouple[i] = indexTouple[i]/dimProd[i];
		for ( int j = i-1 ; j >=0; --j )
			indexTouple[j] -= indexTouple[i]*dimProd[i];
	}
}

template<typename T>
void GridDistribution<T>::general_xyz_to_conseq_column_major(
		std::vector<size_t> const& grid,
		std::vector<size_t> const& indexTouple,
		size_t & conseq_index) const
{
	conseq_index = indexTouple.front();
	for ( int i = 1; i < static_cast<int>(indexTouple.size()) ; ++i)
	{
		conseq_index *= procGridDatak_[i];
		conseq_index += indexTouple[i];
	}
}

template<typename T>
void GridDistribution<T>::general_conseq_to_xyz_column_major(
		std::vector<size_t> const& grid,
		std::vector<size_t> & indexTouple,
		size_t conseq_index) const
{
	//reverse-engineer expressions like
	//	index = (x*dy + y)*dz + z
	//		  = x*dy*dz + y*dz + z
	//for x,y,z using integer division

	int d = static_cast<int>(grid.size());
	std::vector<int> dimProd( d , 1 );
	for (int i = d-2 ; i >= 0 ; --i)
		dimProd[i] = grid[ i+1 ]*dimProd[i+1];

 	indexTouple = std::vector<size_t>( d , conseq_index );
	for ( int i = 0 ;  i < d-1 ; ++i )
	{
		indexTouple[i] = indexTouple[i]/dimProd[i];
		for ( int j = i+1 ;  j < d ; ++j )
			indexTouple[j] -= indexTouple[i]*dimProd[i];
	}
}

template<typename T>
void GridDistribution<T>::redistribute_locally(
		typename auxillary::TemplateTypedefs<T>::scallop_vector & data,
		std::vector<size_t> const& gridIndicesMapOldToNew,
		size_t blockSize) const
{
	if ( buffer_.size() != blockSize )
		buffer_ = std::vector<T>(blockSize);

	//reshuffel the data according to the index map
	std::set<size_t> mappedIndices;
	std::vector< std::vector<std::pair<size_t,size_t> > > copyLoops;
	for ( size_t i = 0; i < gridIndicesMapOldToNew.size(); ++i )
	{
		if ( mappedIndices.find(i) != mappedIndices.end() )
			continue;
		std::vector<std::pair<size_t,size_t> > loop;
		size_t step = i;
		do
		{
			mappedIndices.insert(step);
			if ( step == gridIndicesMapOldToNew[step] )
				break;
			loop.push_back( std::make_pair(step,gridIndicesMapOldToNew[step]) );
			step = gridIndicesMapOldToNew[step];
		} while ( step != i );
		if ( ! loop.empty() )
			loop.pop_back();//we do not need the last copy that closes the loop

		if ( ! loop.empty() )
			copyLoops.push_back( std::move(loop) );
	}
	for ( auto l : copyLoops )
	{
		//fill buffer with the first target
		std::copy( data.begin()+l.front().first*blockSize,
				data.begin()+(l.front().first+1)*blockSize,
				buffer_.begin());

		for ( auto element : l )
		{
			//reshuffel data. The loop makes sure we never overwrite data
			// that is not already copied inside the array or in the buffer.
			std::copy( data.begin()+element.second*blockSize,
					data.begin()+(element.second+1)*blockSize,
					data.begin()+element.first*blockSize);
		}

		//finally, copy the buffer to the last source
		std::copy(buffer_.begin(),
				buffer_.end(),
				data.begin()+l.back().second*blockSize);
	}
}

template<typename T>
std::vector<size_t>
GridDistribution<T>::get_cube_indices_surrounding(
		bool conseqInKGrid,
		std::vector<bT> const& v) const
{
	assert( v.size()== totalGrid_.size() );

	std::vector<size_t> result;
	std::vector<std::pair<size_t,size_t> > minMaxEachDim;
	for ( size_t i = 0 ; i < v.size(); i++ )
	{
		//min, max index in a this->get_nk_dim(i) periodic grid
		bT vfbz = v[i]-std::floor(v[i]);
		int cellstart = static_cast<int>(std::floor(vfbz*totalGrid_[i]));
		int cellend = cellstart + 1;
		cellstart -= (cellstart/totalGrid_[i])*totalGrid_[i];
		cellend -= (cellend/totalGrid_[i])*totalGrid_[i];
		minMaxEachDim.push_back( std::make_pair(static_cast<size_t>(cellstart),static_cast<size_t>(cellend)));
	}

	auto xyz_to_conseq = [&] (std::vector<size_t> const& tuple )
		{
			return conseqInKGrid ? this->k_xyz_to_conseq( tuple )
						: this->R_xyz_to_conseq( tuple );
		};

	//1D
	if ( v.size() == 1 )
	{
		result.push_back( minMaxEachDim[0].first );
		result.push_back( minMaxEachDim[0].second );
		return result;
	}

	//2D counter clock wise starting at min min
	if ( v.size() == 2 )
	{
		// 1 = min, min
		std::vector<size_t> dima = { minMaxEachDim[0].first, minMaxEachDim[1].first };
		result.push_back( xyz_to_conseq( dima ) );

		// 2 = max, min
		dima = { minMaxEachDim[0].second, minMaxEachDim[1].first };
		result.push_back( xyz_to_conseq( dima ) );

		// 3 = max, max
		dima = { minMaxEachDim[0].second, minMaxEachDim[1].second };
		result.push_back( xyz_to_conseq( dima ) );

		// 4 = min, max
		dima = { minMaxEachDim[0].first, minMaxEachDim[1].second };
		result.push_back( xyz_to_conseq( dima ) );
		return result;
	}

	//3D counter clock wise starting at min min bottom, then clock wise top
	if ( v.size() == 3 )
	{
		// 1 = min, min, min
		std::vector<size_t> dima =
			{ minMaxEachDim[0].first, minMaxEachDim[1].first, minMaxEachDim[2].first };
		result.push_back( xyz_to_conseq( dima ) );

		// 2 = max, min , min
		dima = { minMaxEachDim[0].second, minMaxEachDim[1].first, minMaxEachDim[2].first };
		result.push_back( xyz_to_conseq( dima ) );

		// 3 = max, max, min
		dima = { minMaxEachDim[0].second, minMaxEachDim[1].second, minMaxEachDim[2].first };
		result.push_back( xyz_to_conseq( dima ) );

		// 4 = min, max, min
		dima = { minMaxEachDim[0].first, minMaxEachDim[1].second, minMaxEachDim[2].first };
		result.push_back( xyz_to_conseq( dima ) );

		// 5 = min, min, max
		dima = { minMaxEachDim[0].first, minMaxEachDim[1].first, minMaxEachDim[2].second };
		result.push_back( xyz_to_conseq( dima ) );

		// 6 = max, min, max
		dima = { minMaxEachDim[0].second, minMaxEachDim[1].first, minMaxEachDim[2].second };
		result.push_back( xyz_to_conseq( dima ) );

		// 7 = max, max, max
		dima = { minMaxEachDim[0].second, minMaxEachDim[1].second, minMaxEachDim[2].second };
		result.push_back( xyz_to_conseq( dima ) );

		// 8 = min, max, max
		dima = { minMaxEachDim[0].first, minMaxEachDim[1].second, minMaxEachDim[2].second };
		result.push_back( xyz_to_conseq( dima ) );
		return result;
	}

	error_handling::Error(">4D not implement in the k grid!",4);
	return result;
}

} /* namespace parallel */
} /* namespace scallop */
