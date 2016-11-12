/*	This file GeneralizedSusceptibility.hpp is part of scallop.
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
 *  Created on: Nov 11, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/GeneralizedSusceptibility.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
GeneralizedSusceptibility<T>::GeneralizedSusceptibility() :
	MatsubaraImagTimeFourierTransform<T>( /* bool Fermi=*/ false )
{

}

template<typename T>
void GeneralizedSusceptibility<T>::compute_from_gf( GreensFunctionOrbital<T> const& gf )
{

}

} /* namespace gw_flex */
} /* namespace scallop */
