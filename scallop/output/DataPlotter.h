/*	This file DataPlotter.h is part of scallop.
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

#ifndef SCALLOP_OUTPUT_DATAPLOTTER_H_
#define SCALLOP_OUTPUT_DATAPLOTTER_H_

#include "scallop/input/Configuration.h"
#include "scallop/output/TerminalOut.h"
#include "scallop/output/PlotDataProcessing.h"
#include "scallop/output/SingletGapPlotter.h"
#include "scallop/output/SusceptibilityPlotter.h"

namespace scallop
{
namespace output
{

/**
 * This object contains a list of entities to be plotted.
 *
 * It is initialized from the input data and builds the object that in the end
 * plots the data.
 */
class DataPlotter
{
public:

	DataPlotter();

	DataPlotter( input::Configuration const& input );

	template<typename T>
	void plot_per_iteration(output::TerminalOut & msg,
			PlotDataProcessing<T> data);

	template<typename T>
	void plot_end_of_iterations(output::TerminalOut & msg,
			PlotDataProcessing<T> data);

	SusceptibilityPlotter & get_susc_plotter();
private:

	SingletGapPlotter sp_;

	SusceptibilityPlotter spinSuscP_;
};

} /* namespace output */
} /* namespace scallop */

#include "scallop/output/src/DataPlotter.hpp"
#endif /* SCALLOP_OUTPUT_DATAPLOTTER_H_ */
