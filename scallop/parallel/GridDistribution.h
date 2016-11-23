/*	This file GridDistribution.h is part of scallop.
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

#ifndef SCALLOP_PARALLEL_GRIDDISTRIBUTION_H_
#define SCALLOP_PARALLEL_GRIDDISTRIBUTION_H_

#include "scallop/auxillary/TemplateTypedefs.h"
#include <vector>
#include <cstdlib>

namespace scallop
{
namespace parallel
{

/**
 * 	This module keeps track of the grid parallelization in real and reciprocal space.
 *
 *	We implement a parallelization so that each processor handles a given chunk of k vectors.
 *	The distribution is such that for nx grid points in kx direction
 *		| proc 1 | proc 2 | ... | proc N      |
 *	kx =|1|2|..  |nx/N|...| ... |(nx/N)*N-1...|
 *
 *	Then, if we transform to real space, the distribution is along the Ry direction.
 *	This makes it more efficient to synchronize the data.
 */
template<typename T>
class GridDistribution
{
public:

	GridDistribution();

	void distribute_grid( std::vector<size_t> totalGrid );

	/**
	 * Get dimension of space.
	 * @return Dimension of space.
	 */
	size_t get_dim() const;

	/**
	 * Get total space grid.
	 * @return Vector with number of points in each direction x,y,z.
	 */
	std::vector<size_t> const& get_grid() const;

	/**
	 * Get the total number of grid points
	 *
	 * @return total number of grid points
	 */
	size_t get_num_grid() const;

	/**
	 * Get the number of k points this processor handles
	 *
	 * @return Number of k point
	 */
	size_t get_num_k_grid() const;

	std::vector<size_t> const& get_k_grid() const;

	size_t get_num_R_grid() const;

	std::vector<size_t> const& get_R_grid() const;

	size_t get_num_grid_data() const;

	/**
	 * Obtain a tuple of k grid indices for a consecutively ordered index
	 * @param ik	consecutively ordered index
	 * @return		reference to internal buffer with the tuple.
	 */
	std::vector<size_t> & k_conseq_to_xyz( size_t ik ) const;

	/**
	 * Obtain a tuple of R grid indices for a consecutively ordered index
	 * @param iR	consecutively ordered index
	 * @return		reference to internal buffer with the tuple.
	 */
	std::vector<size_t> & R_conseq_to_xyz( size_t iR ) const;

	/**
	 * Convert a processor local consecutive k index to a processor local data grid index.
	 *
	 * @param ik processor local consecutive k index.
	 * @return index w.r.p.t. the data array.
	 */
	size_t k_conseq_local_to_data_conseq( size_t ik ) const;

	/**
	 * Convert a processor local consecutive R index to a processor local data grid index.
	 *
	 * @param iR processor local consecutive R index.
	 * @return index w.r.p.t. the data array.
	 */
	size_t R_conseq_local_to_data_conseq( size_t iR ) const;

	/**
	 * Convert a processor local consecutive k index in to an index tuple in the full grid.
	 *
	 * @param ik processor local consecutive k index.
	 * @return vector of dimension of space containing the indices in the global k grid.
	 */
	std::vector<size_t> & k_conseq_local_to_xyz_total( size_t ik ) const;

	/**
	 * Convert a processor local consecutive R index in to an index tuple in the full grid.
	 *
	 * @param iR processor local consecutive R index.
	 * @return vector of dimension of space containing the indices in the global R grid.
	 */
	std::vector<size_t> & R_conseq_local_to_xyz_total( size_t iR ) const;

	/**
	 * Convert a tuple in the full R grid to a processor local consecutive index.
	 *
	 * In debug mode, this fails if the indices do not correspond to a vector in the
	 * range of this processor.
	 *
	 * @param tuple		xyz indices of the total grid.
	 * @return 			index in the processor local R grid.
	 */
	size_t R_xyz_total_to_conseq_local( std::vector<size_t> tuple ) const;

	/**
	 * Convert a tuple in the full R grid to a consecutive index in the full grid.
	 *
	 * @param tuple		xyz indices of the total R grid.
	 * @return 			index in the full R grid.
	 */
	size_t R_xyz_to_conseq( std::vector<size_t> const& tuple ) const;

	/**
	 * Convert a tuple in the full k grid to a consecutive index in the full grid.
	 *
	 * @param tuple		xyz indices of the total k grid.
	 * @return 			index in the full k grid.
	 */
	size_t k_xyz_to_conseq( std::vector<size_t> const& tuple ) const;

	size_t get_num_grid_total() const;

	size_t get_proc_id_of_gridpoint( size_t igTotal ) const;

	/**
	 * Find the index of the inverse k vector.
	 *
	 * @param ik	index of the vector k
	 * @return		consecutive index of the vector -k (+G with G such that G-k is in the first zone)
	 */
	size_t get_inverse_index_k( size_t ik ) const;

	/**
	 * Find the index of the inverse R vector.
	 *
	 * @param iR	index of the vector R
	 * @return		consecutive index of the vector -R (+T
	 * 				with T such that T-R is in the range of Born-von-Karman periodic boundary conditions)
	 */
	size_t get_inverse_index_R( size_t iR ) const;

	/**
	 * Redistributes data among processors in R/k space to the respective other space.
	 *
	 *	The data ordering in grids is defined in the module as
	 *	(iGx*...)NGPtsy + iGN)
	 *	For both k and R where in k space iGx is processor local while in R space iGy and NGPtsy is processor local.
	 *
	 * @param spaceIsK		When true, the input space of \p data is k,
	 * 						otherwise is the input space of \p data R.
	 * @param data			The container with the data per processor.
	 * 						On input: data to redistributed.
	 * 						On output: redistributed data.
	 * @param blockSize		Number of elements per grid point.
	 */
	void grid_data_transposition(
			bool startFromKSpace,
			typename auxillary::TemplateTypedefs<T>::scallop_vector & data,
			size_t blockSize) const;
private:

	size_t nK_ = 0;

	size_t nR_ = 0;

	size_t nProcGridData_ = 0;

	size_t nGridTotal_ = 0;

	std::vector<size_t> totalGrid_;

	std::vector<size_t> procGridDatak_;

	std::vector<size_t> procGridDataR_;

	std::vector<size_t> procGridk_;

	std::vector<size_t> procGridR_;

	std::vector< std::pair<size_t,size_t> > parallelMapk_;

	std::vector< std::pair<size_t,size_t> > parallelMapR_;

	mutable std::vector<size_t> tupleBuff_;

	mutable std::vector<T> buffer_;

	void distribute_dim( size_t N,
			std::vector< std::pair<size_t,size_t> > & procDimDistribution ) const;

	void general_xyz_to_conseq_row_major(
			std::vector<size_t> const& grid,
			std::vector<size_t> const& indexTouple,
			size_t & conseq_index) const;

	void general_conseq_to_xyz_row_major(
			std::vector<size_t> const& grid,
			std::vector<size_t> & indexTouple,
			size_t conseq_index) const;

	void general_xyz_to_conseq_column_major(
			std::vector<size_t> const& grid,
			std::vector<size_t> const& indexTouple,
			size_t & conseq_index) const;

	void general_conseq_to_xyz_column_major(
			std::vector<size_t> const& grid,
			std::vector<size_t> & indexTouple,
			size_t conseq_index) const;

	void tuple_transpose( std::vector<size_t> & tuple ) const;

	void redistribute_locally(
			typename auxillary::TemplateTypedefs<T>::scallop_vector & data,
			std::vector<size_t> const& gridIndicesMapOldToNew,
			size_t blockSize) const;
};

} /* namespace parallel */
} /* namespace scallop */

#include "scallop/parallel/src/GridDistribution.hpp"
#endif /* SCALLOP_PARALLEL_GRIDDISTRIBUTION_H_ */
