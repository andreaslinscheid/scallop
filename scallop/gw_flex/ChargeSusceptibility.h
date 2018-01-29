/*	This file ChargeSusceptibility.h is part of scallop.
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
 *  Created on: Nov 22, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_CHARGESUSCEPTIBILITY_H_
#define SCALLOP_GW_FLEX_CHARGESUSCEPTIBILITY_H_

#include "scallop/auxillary/TypeMapComplex.h"
#include "scallop/gw_flex/GeneralizedSusceptibility.h"
#include "scallop/output/DataPlotter.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
class ChargeSusceptibility : public GeneralizedSusceptibility<T>
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;

	void compute_from_gf( GreensFunctionOrbital<T> const & GF );

	void copy_charge_part( GeneralizedSusceptibility<T> const& gsust );

	void charge_RPA_enhancement(InteractionMatrix<T> const& interMat);

	void charge_RPA_enhancement(InteractionMatrix<T> const& interMat,
			typename GeneralizedSusceptibility<T>::AdiabaticUpscale & a);

	void charge_effective_interaction(
			InteractionMatrix<T> const& interMat,
			typename GeneralizedSusceptibility<T>::AdiabaticUpscale & a,
			output::DataPlotter & plotter);

	T operator() (size_t ik, size_t iw, size_t m1,  size_t m2) const;

	T & operator() (size_t ik, size_t iw, size_t m1,  size_t m2);

	T operator() (size_t ik, size_t iw, size_t l1, size_t l2, size_t l3, size_t l4) const;

	T & operator() (size_t ik, size_t iw, size_t l1, size_t l2, size_t l3, size_t l4);
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/ChargeSusceptibility.hpp"
#endif /* SCALLOP_GW_FLEX_CHARGESUSCEPTIBILITY_H_ */
