/*	This file PlotDataProcessing.hpp is part of scallop.
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

#include "scallop/output/PlotDataProcessing.h"

namespace scallop
{
namespace output
{

template<typename T>
PlotDataProcessing<T>::PlotDataProcessing( gw_flex::Driver<T> const * driver )
	: driverPtr_(driver)
{

}

template<typename T>
gw_flex::SelfEnergy<T> const &
PlotDataProcessing<T>::get_self_energy() const
{
	return driverPtr_->get_self_energy();
}

template<typename T>
gw_flex::GreensFunctionOrbital<T> const &
PlotDataProcessing<T>::get_greensfunction() const
{
	return driverPtr_->get_greensfunction();
}

} /* namespace output */
} /* namespace scallop */
