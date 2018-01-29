/*	This file SingletGapPlotter.h is part of scallop.
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

#ifndef SCALLOP_OUTPUT_SINGLETGAPPLOTTER_H_
#define SCALLOP_OUTPUT_SINGLETGAPPLOTTER_H_

#include "scallop/output/PlotDataProcessing.h"
#include "scallop/output/LowestMatsFreqOnSpaceGridPlotter.h"

namespace scallop
{
namespace output
{

class SingletGapPlotter : public LowestMatsFreqOnSpaceGridPlotter
{
public:

	template<typename T>
	void plot(PlotDataProcessing<T> data);
};

} /* namespace output */
} /* namespace scallop */

#include "scallop/output/src/SingletGapPlotter.hpp"
#endif /* SCALLOP_OUTPUT_SINGLETGAPPLOTTER_H_ */
