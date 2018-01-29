/*	This file FormatDescriptor.h is part of scallop.
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

#ifndef SCALLOP_OUTPUT_FORMATDESCRIPTOR_H_
#define SCALLOP_OUTPUT_FORMATDESCRIPTOR_H_

#include <string>
#include <vector>
#include <complex>

namespace scallop
{
namespace output
{

class FormatDescriptor
{
public:

	typedef enum
	{
		defaultPlotProg = 0,
		ParaView = 1,
		Gnuplot = 2,
	} PlotProgram;

	void setup(
			std::string filename,
			std::string quantityName,
			PlotProgram pp = defaultPlotProg);

	PlotProgram get_plot_program() const;

	std::string get_quantity_label() const;

	bool do_plot_per_iteration() const;

	std::string get_filename() const;

	void set_active_plotprogram(PlotProgram pp);
private:

	PlotProgram activePP_ = PlotProgram::Gnuplot;

	std::string fileName_;

	std::string quantityName_;
};

} /* namespace output */
} /* namespace scallop */

#endif /* SCALLOP_OUTPUT_FORMATDESCRIPTOR_H_ */
