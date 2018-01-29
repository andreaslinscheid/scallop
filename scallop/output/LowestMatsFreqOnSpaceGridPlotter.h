/*	This file LowestMatsFreqOnSpaceGridPlotter.h is part of scallop.
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
 *  Created on: Jan 21, 2017
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_OUTPUT_LOWESTMATSFREQONSPACEGRIDPLOTTER_H_
#define SCALLOP_OUTPUT_LOWESTMATSFREQONSPACEGRIDPLOTTER_H_

#include "scallop/gw_flex/MatsubaraImagTimeFourierTransform.h"
#include "scallop/output/FormatDescriptor.h"
#include <string>

namespace scallop
{
namespace output
{

/**
 * Designed to handle the output to disk of objects to be plotted.
 */
class LowestMatsFreqOnSpaceGridPlotter : public FormatDescriptor
{
public:

	void plot_2D_grid_data(
			std::vector<float> data,
			std::vector<float> gridInY,
			std::vector<float> gridInX) const;

};

} /* namespace output */
} /* namespace scallop */

#include "scallop/output/src/LowestMatsFreqOnSpaceGridPlotter.hpp"
#endif /* SCALLOP_OUTPUT_LOWESTMATSFREQONSPACEGRIDPLOTTER_H_ */
