/*	This file SpinSusceptibility.h is part of scallop.
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

#ifndef SCALLOP_GW_FLEX_SPINSUSCEPTIBILITY_H_
#define SCALLOP_GW_FLEX_SPINSUSCEPTIBILITY_H_

#include "scallop/gw_flex/GeneralizedSusceptibility.h"
#include "scallop/output/DataPlotter.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
class SpinSusceptibility : public GeneralizedSusceptibility<T>
{
public:

	void compute_from_gf( GreensFunctionOrbital<T> const & GF);

	using GeneralizedSusceptibility<T>::operator();

	void spin_RPA_enhancement(
			InteractionMatrix<T> const& interMat);

	void spin_RPA_enhancement(
			InteractionMatrix<T> const& interMat,
			typename GeneralizedSusceptibility<T>::AdiabaticUpscale & a);

	void spin_effective_interaction(
			InteractionMatrix<T> const& interMat,
			typename GeneralizedSusceptibility<T>::AdiabaticUpscale & a,
			output::DataPlotter & plotter);
private:

	size_t nChannels_ = 3;
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/SpinSusceptibility.hpp"
#endif /* SCALLOP_GW_FLEX_SPINSUSCEPTIBILITY_H_ */
