/*	This file GnuplotBinaryData2DFile.h is part of scallop.
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

#ifndef SCALLOP_OUTPUT_GNUPLOTBINARYDATA2DFILE_H_
#define SCALLOP_OUTPUT_GNUPLOTBINARYDATA2DFILE_H_

#include <vector>
#include <string>

namespace scallop
{
namespace output
{

class GnuplotBinaryData2DFile
{
public:

	/**
	 * Write a binary file in parallel with a format suitable for gnuplot.
	 *
	 * The data must be distributed such that proc 1 has all y for the smallest x
	 * then proc 2 has all y for next smallest x and so on.
	 *
	 * @param filename	The name of the file to be written.
	 * @param data		The data written into the file in the processor order as explained above.
	 * @param gridInY	The grid points in y
	 * @param gridInX	The grid points in x
	 */
	void write_gnuplot_binary_2d_data_file(
			std::string const& filename,
			std::vector<float> data,
			std::vector<float> gridInY,
			std::vector<float> gridInX) const;
};

} /* namespace output */
} /* namespace scallop */

#endif /* SCALLOP_OUTPUT_GNUPLOTBINARYDATA2DFILE_H_ */
