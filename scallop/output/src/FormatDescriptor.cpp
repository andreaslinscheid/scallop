/*	This file FormatDescriptor.cpp is part of scallop.
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

#include "scallop/output/FormatDescriptor.h"

namespace scallop
{
namespace output
{

void FormatDescriptor::setup(
		std::string filename,
		std::string quantityName,
		PlotProgram pp)
{
	fileName_ = std::move(filename);
	quantityName_ = std::move(quantityName);
	activePP_ = pp;
}

bool FormatDescriptor::do_plot_per_iteration() const
{
	//todo implement another switch to allow to only to plot in the end
	return (not fileName_.empty() );
}

std::string FormatDescriptor::get_filename() const
{
	return fileName_;
}

void FormatDescriptor::set_active_plotprogram(PlotProgram pp)
{
	activePP_ = pp;
}

FormatDescriptor::PlotProgram
FormatDescriptor::get_plot_program() const
{
	return activePP_;
}

std::string FormatDescriptor::get_quantity_label() const
{
	return quantityName_;
}

} /* namespace output */
} /* namespace scallop */
