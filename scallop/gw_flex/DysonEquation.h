/*	This file DysonEquation.h is part of scallop.
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
 *  Created on: Dec 16, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_DYSONEQUATION_H_
#define SCALLOP_GW_FLEX_DYSONEQUATION_H_

#include "scallop/gw_flex/GreensFunctionOrbital.h"
#include "scallop/auxillary/TemplateTypedefs.h"
#include "scallop/gw_flex/KohnShamBandStructure.h"
#include "scallop/gw_flex/SelfEnergy.h"

namespace scallop
{
namespace gw_flex
{

class DysonEquation
{
public:

	template<typename T>
	void solve_by_inversion(
			GreensFunctionOrbital<T> & g,
			typename auxillary::TemplateTypedefs<T>::scallop_vector const& MatsFreqs,
			KohnShamBandStructure<T> const& elstr,
			SelfEnergy<T> const& sigma) const;
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/DysonEquation.hpp"
#endif /* SCALLOP_GW_FLEX_DYSONEQUATION_H_ */
