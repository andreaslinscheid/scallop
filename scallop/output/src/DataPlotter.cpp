/*	This file DataPlotter.cpp is part of scallop.
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
#include "scallop/output/SingletGapPlotter.h"

namespace scallop
{
namespace output
{

DataPlotter::DataPlotter()
{

}

DataPlotter::DataPlotter( input::Configuration const& input )
{
	//Check if we plot the Green's function's spectral function.
	if ( not input.get_f_spec().empty() )
	{

	}

	//Check if we plot the singlet gap as a function of the k grid.
	if ( not input.get_f_gap_s().empty() )
	{
		sp_.setup( input.get_f_gap_s(), "${\\Delta}^{s}_{\\bf{k}}$" );
	}

	//Check if we plot the spin susceptibility
	if ( not input.get_f_spin_sust_stat().empty() )
	{
		spinSuscP_.configure_spin(true);
		spinSuscP_.setup( input.get_f_spin_sust_stat(), "${\\rm{max EV}(\\chi}^{s,\\rm{RPA}}_{\\bf{k}})$" );
	}
}

SusceptibilityPlotter &
DataPlotter::get_susc_plotter()
{
	return spinSuscP_;
}

} /* namespace output */
} /* namespace scallop */
