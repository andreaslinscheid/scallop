/*	This file SpinSusceptibility.hpp is part of scallop.
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
 *  Created on: Dec 15, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/SpinSusceptibility.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
void SpinSusceptibility<T>::compute_from_gf( GreensFunctionOrbital<T> const & GF)
{
	GeneralizedSusceptibility<T>::compute_from_gf(GF,nChannels_,/*offset=*/1);
}

template<typename T>
void SpinSusceptibility<T>::spin_RPA_enhancement(
		InteractionMatrix<T> const& interMat)
{
	GeneralizedSusceptibility<T>::RPA_enhancement(interMat);
}

template<typename T>
void SpinSusceptibility<T>::spin_effective_interaction(
		InteractionMatrix<T> const& interMat,
		typename GeneralizedSusceptibility<T>::AdiabaticUpscale & a,
		output::DataPlotter & plotter)
{
	if ( plotter.get_susc_plotter().do_plot_static() )
	{
		output::TerminalOut msg;
		msg << "Plotting largest eigenvalue of the static spin susceptibility ...";
	}
	GeneralizedSusceptibility<T>::RPA_enhancement(interMat,a,plotter);
}

} /* namespace gw_flex */
} /* namespace scallop */
