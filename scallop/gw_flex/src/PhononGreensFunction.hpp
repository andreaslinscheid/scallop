/*	This file PhononGreensFunction.hpp is part of scallop.
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

#include "scallop/gw_flex/PhononGreensFunction.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
PhononGreensFunction<T>::PhononGreensFunction() :
		MatsubaraImagTimeFourierTransform<T>( /* fermionic = */ false )
{

}

template<typename T>
void PhononGreensFunction<T>::initialize(
		size_t dimImTime,
		std::vector<size_t> gridDims,
		size_t numModes,
		bool initialInTimeDomain,
		bool initialInReciprocalDomain,
		typename auxillary::TemplateTypedefs<T>::scallop_vector const& data)
{

	this->initialize_layout_phonon_prop(numModes);

	MatsubaraImagTimeFourierTransform<T>::initialize(
			dimImTime,
			gridDims,
			numModes*numModes,
			initialInTimeDomain,
			initialInReciprocalDomain,
			data);
}

template<typename T>
T PhononGreensFunction<T>::operator() (size_t iq, size_t iw, size_t nu, size_t mu) const
{
	return *(this->read_phs_grid_ptr_block(iq,iw) + memory_layout_phonon_prop(nu,mu));
}

template<typename T>
T & PhononGreensFunction<T>::operator() (size_t iq, size_t iw, size_t nu, size_t mu)
{
	return *(this->write_phs_grid_ptr_block(iq,iw) + memory_layout_phonon_prop(nu,mu));
}

} /* namespace gw_flex */
} /* namespace scallop */
