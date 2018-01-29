/*	This file PlotDataProcessing.h is part of scallop.
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

#ifndef SCALLOP_OUTPUT_PLOTDATAPROCESSING_H_
#define SCALLOP_OUTPUT_PLOTDATAPROCESSING_H_

#include "scallop/gw_flex/SelfEnergy.h"
#include "scallop/gw_flex/GreensFunctionOrbital.h"

namespace scallop
{
namespace gw_flex
{
template<typename T>
class Driver; //forward declare
};
};

namespace scallop
{
namespace output
{

template<typename T>
class PlotDataProcessing
{
public:

	PlotDataProcessing( gw_flex::Driver<T> const * driver );

	gw_flex::SelfEnergy<T> const & get_self_energy() const;

	gw_flex::GreensFunctionOrbital<T> const & get_greensfunction() const;

private:

	gw_flex::Driver<T> const * driverPtr_;
};

} /* namespace output */
} /* namespace scallop */

#include "scallop/output/src/PlotDataProcessing.hpp"
#endif /* SCALLOP_OUTPUT_PLOTDATAPROCESSING_H_ */
