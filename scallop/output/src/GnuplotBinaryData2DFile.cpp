/*	This file GnuplotBinaryData2DFile.cpp is part of scallop.
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
 *  Created on: Feb 6, 2017
 *      Author: A. Linscheid
 */

#include "scallop/output/GnuplotBinaryData2DFile.h"
#include "scallop/parallel/MPIModule.h"
#include <assert.h>

namespace scallop
{
namespace output
{

void GnuplotBinaryData2DFile::write_gnuplot_binary_2d_data_file(
		std::string const& filename,
		std::vector<float> data,
		std::vector<float> gridInY,
		std::vector<float> gridInX) const
{
	assert( gridInX.size()*gridInY.size() == data.size() );

	//data has to be padded with the grid so that gnuplot can plot it
	//this means on the first processor we need to add space for the y grid.
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	std::vector< float > dataForGP( data.size()+gridInX.size()
							+ (mpi.get_mpi_me()==0?1+gridInY.size():0));

	size_t pos = 0;
	if ( mpi.get_mpi_me() == 0 )
	{
		dataForGP[pos++] = gridInY.size();

		for(auto gy:gridInY)
			dataForGP[pos++] = gy;
	}

	for ( size_t ix = 0; ix < gridInX.size(); ++ix)
	{
		dataForGP[pos++] = gridInX[ix];

		for(size_t i = 0; i < gridInY.size(); ++i)
			dataForGP[pos++] = data[ix*gridInY.size()+i];
	}

	mpi.save_file_parallel( dataForGP, filename);
}

} /* namespace output */
} /* namespace scallop */
