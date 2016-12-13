/*	This file IrregularGridDistribution.h is part of scallop.
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

#ifndef SCALLOP_PARALLEL_IRREGULARGRIDDISTRIBUTION_H_
#define SCALLOP_PARALLEL_IRREGULARGRIDDISTRIBUTION_H_

#include "scallop/parallel/GridDistribution.h"

namespace scallop
{
namespace parallel
{

template<typename T>
class IrregularGridDistribution
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;
	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector V;
	typedef typename auxillary::TemplateTypedefs<bT>::scallop_vector VbT;

	void distribute_pts(
			bool inKSpace,
			VbT points,
			GridDistribution<T> regularGrid );

	std::vector<bT> get_vector( size_t localIndex ) const;

	size_t proc_local_to_total_index(size_t ig) const;

	size_t get_n_pts_total() const;

	size_t get_n_pts() const;

	size_t get_dim() const;

	/**
	 * Obtain a list of all corner points of cubes in the regular grid that enclose the current set of points.
	 *
	 * @param list	On return list contains all consecutive grid indices in the full grid that are need for the current set of points
	 */
	void get_list_required_regular_grid_indices( std::vector< size_t > & list) const;

	/**
	 * Communicate data on the regular grid used such that interpolation of a processor local operation.
	 *
	 * @param dRegularWedge	The data on the regular grid used for distribute_pts, distributed according to GridDistribution
	 * @param dReqGrid		Per processor, the data at each grid point specified by
	 * 						::get_list_required_regular_grid_indices(std::vector< size_t > & list) const
	 */
	template<class DI, class DO>
	void proc_sync_data(
			bool isInKSpace,
			DI const& dRegularWedge,
			DO & dReqGrid,
			size_t blocksizePerGridPt) const;

	template<class DI, class DO>
	void linear_interpolate_data(
			DI const& data,
			DO & interpolated_data,
			size_t blocksizePerGridPt) const;

	template<typename FI,typename FO>
	void linear_interpolate_data(
			FI const * data,
			FO * interpolated_data,
			size_t blocksizePerGridPt) const;

	GridDistribution<T> const& get_spaceGrid_proc() const;
private:

	struct GridCube
	{
		GridCube( std::vector<size_t> cornerPoints )
				: cornerPoints_( std::move(cornerPoints) ){};

		std::vector<size_t> cornerPoints_;
		mutable std::vector<size_t> containedIrregularPts_;

		bool operator< (GridCube const& other) const;
	};

	size_t nPts_ = 0;

	size_t nPtsTotal_ = 0;

	size_t dim_ = 0;

	VbT ptsCoordinates_;

	GridDistribution<T> gridDistr_;

	std::vector< size_t > conseqPtsRegularGrid_;

	std::vector< GridCube > cubes_;

	std::vector< size_t > occurConseqRegularGridIndices_;

	std::vector<int> procOffset_;

	template<typename F>
	F bilinear_interpol(bT x, bT y, F f00, F f10, F f11, F f01) const;

	template<typename F>
	F trilinear_interpol(bT x, bT y, bT z, F f000, F f100, F f110, F f010, F f001, F f101, F f111, F f011) const;
};

} /* namespace parallel */
} /* namespace scallop */

#include "scallop/parallel/src/IrregularGridDistribution.hpp"
#endif /* SCALLOP_PARALLEL_IRREGULARGRIDDISTRIBUTION_H_ */
