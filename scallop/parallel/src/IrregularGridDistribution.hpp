/*	This file IrregularGridDistribution.hpp is part of scallop.
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
 *  Created on: Nov 28, 2016
 *      Author: A. Linscheid
 */

#include "scallop/parallel/IrregularGridDistribution.h"
#include "scallop/parallel/MPIModule.h"
#include "scallop/auxillary/DataRearrangement.h"
#include <map>
#include <set>

namespace scallop
{
namespace parallel
{

template<typename T>
void IrregularGridDistribution<T>::distribute_pts(
		bool inKSpace,
		VbT points,
		GridDistribution<T> regularGrid )
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	gridDistr_ = std::move(regularGrid);
	dim_ = gridDistr_.get_grid().size();
	std::vector<std::pair<size_t,size_t> > distrb;
	nPtsTotal_ = points.size()/dim_;
	gridDistr_.distribute_dim(nPtsTotal_,distrb);
	int first = static_cast<int>(distrb[ mpi.get_mpi_me() ].first);
	int last =  static_cast<int>(distrb[ mpi.get_mpi_me() ].second);
	procOffset_ = std::vector<int>( mpi.get_nproc() , 0 );
	for ( size_t ip = 1 ; ip < mpi.get_nproc(); ++ip)
		procOffset_[ip] = static_cast<int>(distrb[ip-1].second);

	nPts_ = static_cast<size_t>( last - first );
	ptsCoordinates_ = VbT(nPts_*dim_);
	//We keep the coordinates strictly in the zone [0,1[
	for ( size_t i = 0 ; i <nPts_*dim_; ++i )
		ptsCoordinates_[i] = points[first*dim_+i]-std::floor(points[first*dim_+i]);

	//here we compute all the indices of corner points of cubes in the regular grid
	// and accumulate all irregular grid points in their respective cube
	// and save the cubes for later use
	std::set< GridCube > cubeSet;
	std::vector<bT> v(dim_);
	auto lastUsed = cubeSet.end();
	for ( size_t ip = 0 ; ip < nPts_; ++ip)
	{
		auto it = ptsCoordinates_.begin()+ip*dim_;
		std::copy(it,it+dim_,v.begin() );
		GridCube c( (gridDistr_.get_cube_indices_surrounding(inKSpace,v)) );
		// we hint the previous location for insert since we expect the irregular grid points
		//	to be close together
		auto ret = cubeSet.insert(lastUsed,c);
		ret->containedIrregularPts_.push_back(ip);
		lastUsed = ret;
	}
	cubes_ = std::vector< GridCube >( cubeSet.begin(), cubeSet.end() );

	//Now we determine all appearing regular grid indices
	std::set<size_t> conseqRegularGridIndices;
	for ( auto cube : cubes_)
		conseqRegularGridIndices.insert( cube.cornerPoints_.begin(), cube.cornerPoints_.end() );

	//and here we cast them into an internal mapping
	// If the code below is complete, cubes_ contain per
	//	grid point, indices i \in [0...#appearing grid indices]
	//	and conseqPtsRegularGrid_ that resolves i to the consecutive ordered index.
	conseqPtsRegularGrid_ = std::vector<size_t>( conseqRegularGridIndices.size() );
	std::map<size_t,size_t> invConseqPtsRegularGrid;
	size_t i = 0;
	for ( auto ig : conseqRegularGridIndices )
	{
		invConseqPtsRegularGrid[ig] = i; // guaranteed to be unique by set
		conseqPtsRegularGrid_[i++]= ig;
	}
	size_t nPtsThisProc = 0;
	for ( auto & cube : cubes_)
	{
		for ( auto & coi : cube.cornerPoints_ )
		{
			auto it = invConseqPtsRegularGrid.find( coi );
			if ( it == invConseqPtsRegularGrid.end() )
				error_handling::Error("Internal programming error!");
			coi = it->second;
		}
		nPtsThisProc += cube.containedIrregularPts_.size();
	}

	if ( this->get_n_pts() != nPtsThisProc)
		error_handling::Error("We have lost some points. Internal programming error!");
}

template<typename T>
size_t IrregularGridDistribution<T>::get_n_pts_total() const
{
	return nPtsTotal_;
}

template<typename T>
size_t IrregularGridDistribution<T>::proc_local_to_total_index(size_t ig) const
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	return static_cast<size_t>( static_cast<int>(ig) + procOffset_[ mpi.get_mpi_me() ] );
}

template<typename T>
size_t IrregularGridDistribution<T>::get_n_pts() const
{
	return nPts_;
}

template<typename T>
void IrregularGridDistribution<T>::get_list_required_regular_grid_indices(
		std::vector< size_t > & list) const
{
	list = conseqPtsRegularGrid_;
}

template<typename T>
template<typename FI,typename FO>
void IrregularGridDistribution<T>::linear_interpolate_data(
	FI const * data,
	FO * interpolated_data,
	size_t blocksizePerGridPt) const
{
	size_t nB = blocksizePerGridPt;

	auto k = std::vector<bT>(dim_);
	auto fill_k = [&] (size_t irredKIndex)
	{
		assert( irredKIndex < ptsCoordinates_.size()/dim_ );
		auto it = ptsCoordinates_.begin()+irredKIndex*dim_;
		std::copy( it, it+dim_, k.begin() );
	};

	typedef typename auxillary::TypeMapComplex<FI>::type real_FI;

	std::vector<bT> kcell1, kcell2;
	for ( size_t icube = 0; icube < cubes_.size(); ++icube )
	{
		//use one representative k vector to compute the cell vectors
		fill_k( cubes_[icube].containedIrregularPts_.front() );
		gridDistr_.get_cell_vectors(k,kcell1, kcell2);

		assert( cubes_[icube].cornerPoints_.size() == std::pow(2,dim_));
		//1D
		if ( dim_ == 1 )
		{
			auto g1it = &( data[ nB*cubes_[icube].cornerPoints_[0]] );
			auto g2it = &( data[ nB*cubes_[icube].cornerPoints_[1]] );

			for ( size_t ikirred : cubes_[icube].containedIrregularPts_ )
			{
				fill_k(ikirred);
				bT x = (k[0] - kcell1[0])/(kcell2[0] - kcell1[0]);
				for ( size_t ib= 0 ; ib < nB; ++ib)
				{
					interpolated_data[ikirred*nB+ib]= FO(g1it[ib] * real_FI(1. - x) + g2it[ib] * real_FI(x));
					assert( interpolated_data[ikirred*nB+ib] == interpolated_data[ikirred*nB+ib]);
				}
			}
		}
		//2D
		if ( dim_ == 2 )
		{
			auto g1it = &( data[ nB*cubes_[icube].cornerPoints_[0]] );
			auto g2it = &( data[ nB*cubes_[icube].cornerPoints_[1]] );
			auto g3it = &( data[ nB*cubes_[icube].cornerPoints_[2]] );
			auto g4it = &( data[ nB*cubes_[icube].cornerPoints_[3]] );

			for ( size_t ikirred : cubes_[icube].containedIrregularPts_ )
			{
				fill_k(ikirred);
				bT x = (k[0] - kcell1[0])/(kcell2[0] - kcell1[0]);
				bT y = (k[1] - kcell1[1])/(kcell2[1] - kcell1[1]);
				for ( size_t ib= 0 ; ib < nB; ++ib)
				{
					interpolated_data[ikirred*nB+ib]=  FO(this->bilinear_interpol( x , y,
								g1it[ib], g2it[ib],
								g3it[ib], g4it[ib] ));
					assert( interpolated_data[ikirred*nB+ib] == interpolated_data[ikirred*nB+ib]);
				}
			}
		}
		//3D
		if ( dim_ == 3 )
		{
			auto g1it = &( data[ nB*cubes_[icube].cornerPoints_[0]] );
			auto g2it = &( data[ nB*cubes_[icube].cornerPoints_[1]] );
			auto g3it = &( data[ nB*cubes_[icube].cornerPoints_[2]] );
			auto g4it = &( data[ nB*cubes_[icube].cornerPoints_[3]] );
			auto g5it = &( data[ nB*cubes_[icube].cornerPoints_[4]] );
			auto g6it = &( data[ nB*cubes_[icube].cornerPoints_[5]] );
			auto g7it = &( data[ nB*cubes_[icube].cornerPoints_[6]] );
			auto g8it = &( data[ nB*cubes_[icube].cornerPoints_[7]] );

			for ( size_t ikirred : cubes_[icube].containedIrregularPts_ )
			{
				fill_k(ikirred);
				bT x = (k[0] - kcell1[0])/(kcell2[0] - kcell1[0]);
				bT y = (k[1] - kcell1[1])/(kcell2[1] - kcell1[1]);
				bT z = (k[2] - kcell1[2])/(kcell2[2] - kcell1[2]);
				for ( size_t ib= 0 ; ib < nB; ++ib)
				{
					interpolated_data[ikirred*nB+ib] = FO(this->trilinear_interpol( x , y, z,
							g1it[ib], g2it[ib],
							g3it[ib], g4it[ib],
							g5it[ib], g6it[ib],
							g7it[ib], g8it[ib] )) ;
					assert( interpolated_data[ikirred*nB+ib] == interpolated_data[ikirred*nB+ib]);
				}
			}
		}
	}
}

template<typename T>
template<class DI, class DO>
void IrregularGridDistribution<T>::linear_interpolate_data(
	DI const& data,
	DO & interpolated_data,
	size_t blocksizePerGridPt) const
{
	assert( data.size() == blocksizePerGridPt*conseqPtsRegularGrid_.size() );

	if ( interpolated_data.empty() )
		interpolated_data = DO( blocksizePerGridPt*this->get_n_pts() );

	this->linear_interpolate_data(data.data(),interpolated_data.data(),blocksizePerGridPt);
}


template<typename T>
bool IrregularGridDistribution<T>::GridCube::operator< (GridCube const& other) const
{
	for ( size_t i = 0 ; i < cornerPoints_.size(); ++i)
	{
		if ( cornerPoints_[i] < other.cornerPoints_[i] )
			return true;
		if ( other.cornerPoints_[i] < cornerPoints_[i] )
			return false;
	}
	return false;
};

template<typename T>
template<typename F>
F IrregularGridDistribution<T>::bilinear_interpol(
		bT x, bT y,
		F f00, F f10, F f11, F f01) const
{
	typedef typename auxillary::TypeMapComplex<F>::type real_t;
	F a = f00 * real_t(1. - x) + f10 * real_t(x);
	F b = f01 * real_t(1. - x) + f11 * real_t(x);
	return a * real_t(1. - y) + b * real_t(y);
}

template<typename T>
template<typename F>
F IrregularGridDistribution<T>::trilinear_interpol(
		bT x, bT y, bT z,
		F f000, F f100, F f110, F f010, F f001, F f101, F f111, F f011) const
{
	typedef typename auxillary::TypeMapComplex<F>::type real_t;
	F e = this->bilinear_interpol(x, y, f000, f100, f110, f010);
	F f = this->bilinear_interpol(x, y, f001, f101, f111, f011);
	return e * real_t( 1. - z) + f * real_t(z);
}

template<typename T>
template<class DI, class DO>
void IrregularGridDistribution<T>::proc_sync_data(
		bool isInKSpace,
		DI const& dRegularWedge,
		DO & dReqGrid,
		size_t blocksizePerGridPt) const
{
	// The algorithm is that every processor computes a list of required total grid points.
	// We send this data to every processor so that at each processor we know which data will
	// be used by somebody.
	//	Then, we determine which processor has these grid points and build a
	//			0. send buffer with the local data for all required grid points.
	//			1. list with number of elements to be send to proc 0,1...N
	//	We communicate the list with elements, so that we know how much data is coming from where
	//	so that we can allocate a receive buffer of appropriate size and displacements.
	//	We communicate a list for each processor that tells for each point what its index was
	//	We rearrange the data such that its ordering obeys the result of get_list_required_regular_grid_indices
	MPIModule const& mpi = parallel::MPIModule::get_instance();
	size_t nP = mpi.get_nproc();

	std::vector<size_t> list,requestlistForAll;
	this->get_list_required_regular_grid_indices(list);
	mpi.all_gather( requestlistForAll, list );
	std::vector<size_t> requestNum( nP, 0 );
	std::vector<size_t> sendcount( nP, 0 );
	std::vector< std::vector<size_t> > indexMapReceivGridIndex( nP );
	requestNum[ mpi.get_mpi_me() ] = list.size();
	mpi.sum(requestNum);
	size_t c = 0;
	for ( size_t iproc = 0 ; iproc < nP; ++iproc )
	{
		for ( size_t id = 0 ; id < requestNum[iproc]; ++id)
		{
			size_t originProc = gridDistr_.get_proc_index( isInKSpace, requestlistForAll[c] );
			if ( originProc == mpi.get_mpi_me() )
			{
				sendcount[iproc]++;
			}
			c++;
		}
	}

#ifdef DEBUG_BUILD
	//Check that the total number of send and received items agrees
	size_t n = 0, np = 0;
	for ( size_t iproc = 0 ; iproc < nP; ++iproc )
		n += sendcount[iproc];

	np += list.size();
	mpi.sum(n);
	mpi.sum(np);
	if ( ! (n == np) )
		error_handling::Error( std::string("Internal error: received and send elements need to sum up to the same number but they are")
				+std::to_string(n)+std::string(" and ")+std::to_string(np));

	size_t nrgrid = ( isInKSpace ? gridDistr_.get_num_k_grid():gridDistr_.get_num_R_grid());
	assert( dRegularWedge.size() == blocksizePerGridPt*nrgrid );
#endif

	std::vector<size_t> displ(nP,0);
	for ( int iproc = 1 ; iproc < static_cast<int>(nP); ++iproc )
		displ[iproc] = displ[iproc-1]+sendcount[iproc-1];

	//Fill the send buffer and the indexMap that will be shipped with the data
	size_t nElemSend =  displ.back() + sendcount.back();
	std::vector<size_t> indexMap( nElemSend );
	std::vector<size_t> incr(nP,0);
	typename std::remove_reference<decltype(dReqGrid)>::type sendbuffer( nElemSend*blocksizePerGridPt );
	c = 0;
	for ( size_t ipDest = 0 ; ipDest < nP; ++ipDest )
	{
		for ( size_t id = 0 ; id < requestNum[ipDest]; ++id)
		{
			if ( gridDistr_.get_proc_index( isInKSpace, requestlistForAll[c] )
					== mpi.get_mpi_me() )
			{
				size_t buffIndex = displ[ipDest]+incr[ipDest];
				assert( incr[ipDest] < sendcount[ ipDest ] );
				assert( buffIndex < nElemSend );

				indexMap[ buffIndex ] = requestlistForAll[c];

				size_t procLocalGridIndex = gridDistr_.conseq_full_to_conseq_local( isInKSpace, requestlistForAll[c] );
				assert( procLocalGridIndex*blocksizePerGridPt < dRegularWedge.size() );

				auto it = sendbuffer.begin()+ buffIndex*blocksizePerGridPt;
				auto itSource = dRegularWedge.begin()+ procLocalGridIndex*blocksizePerGridPt;

				std::copy( itSource, itSource+blocksizePerGridPt, it);
				incr[ipDest]++;
			}
			c++;
		}
	}

	auto sendcountBlock =sendcount;
	for ( auto & s : sendcountBlock)
		s *= blocksizePerGridPt;

#ifndef NDEBUG
	for ( auto s: sendbuffer)
		assert( s == s);
#endif

	mpi.all_to_allv( sendbuffer, dReqGrid, sendcountBlock );

	std::vector<size_t> recvIndexMap;
	mpi.all_to_allv( indexMap, recvIndexMap, sendcount );

	std::map<size_t,size_t> oldToNewMap;
	for ( size_t iOld = 0 ; iOld < recvIndexMap.size(); ++iOld)
		oldToNewMap[recvIndexMap[iOld]] = iOld;

	auto it =oldToNewMap.begin();
	std::vector<size_t> indexMapOldToNewArrangement(oldToNewMap.size());
	for ( auto &e : indexMapOldToNewArrangement )
	{
		e = it->second;
		++it;
	}

	auxillary::DataRearrangement< typename std::remove_reference<decltype(dReqGrid)>::type > redistr;
	redistr.redistribute_locally( dReqGrid, indexMapOldToNewArrangement, blocksizePerGridPt );
}

template<typename T>
std::vector< typename IrregularGridDistribution<T>::bT >
IrregularGridDistribution<T>::get_vector( size_t localIndex ) const
{
	auto it = ptsCoordinates_.begin()+localIndex*dim_;
	std::vector<bT> result(it,it+dim_);
	return result;
}

template<typename T>
GridDistribution<T> const&
IrregularGridDistribution<T>::get_spaceGrid_proc() const
{
	return gridDistr_;
}

} /* namespace parallel */
} /* namespace scallop */
