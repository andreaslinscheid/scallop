/*	This file PhononGreensFunction.h is part of scallop.
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
 *  Created on: Nov 14, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_PHONONGREENSFUNCTION_H_
#define SCALLOP_GW_FLEX_PHONONGREENSFUNCTION_H_

#include "scallop/gw_flex/MatsubaraImagTimeFourierTransform.h"
#include "scallop/gw_flex/MemoryLayout.h"


namespace scallop
{
namespace gw_flex
{

template<typename T>
class PhononGreensFunction : public MatsubaraImagTimeFourierTransform<T>, private MemoryLayout
{
public:

	typedef typename auxillary::TypeMapComplex<T>::type bT;

	PhononGreensFunction();

	void initialize(
			size_t dimImTime,
			std::vector<size_t> gridDims,
			size_t numModes,
			bool initialInTimeDomain,
			bool initialInReciprocalDomain,
			typename auxillary::TemplateTypedefs<T>::scallop_vector const& data);

	T operator() (size_t iq, size_t iw, size_t nu, size_t mu) const;

	T & operator() (size_t iq, size_t iw, size_t nu, size_t mu);

};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/PhononGreensFunction.hpp"
#endif /* SCALLOP_GW_FLEX_PHONONGREENSFUNCTION_H_ */
