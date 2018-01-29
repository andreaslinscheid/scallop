/*	This file DataPlotter.hpp is part of scallop.
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

#include "scallop/output/DataPlotter.h"

namespace scallop
{
namespace output
{

template<typename T>
void DataPlotter::plot_per_iteration(output::TerminalOut & msg,
		PlotDataProcessing<T> data)
{
	if ( sp_.do_plot_per_iteration() )
	{
		msg << "\tPlotting singlet gap ...";
		sp_.plot(data);
	}
}

template<typename T>
void DataPlotter::plot_end_of_iterations(output::TerminalOut & msg,
		PlotDataProcessing<T> data)
{

}

} /* namespace output */
} /* namespace scallop */
