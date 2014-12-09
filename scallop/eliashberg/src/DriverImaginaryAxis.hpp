/*	This file DriverImaginaryAxis.cpp is part of scallop.
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
 *  Created on: Nov 25, 2014
 *      Author: Andreas Linscheid
 */

#include "scallop/eliashberg/DriverImaginaryAxis.h"

namespace scallop {
namespace eliashberg {

template<typename T>
DriverImaginaryAxis<T>::DriverImaginaryAxis() :
		_diagonalCouplingUpSpin(),
		_diagonalCouplingDownSpin(),
		_offDiagonalCoupling(),
		_singleRunDriver(_diagonalCouplingUpSpin,_diagonalCouplingDownSpin,_offDiagonalCoupling) { };

template<typename T>
void DriverImaginaryAxis<T>::iterate() {
	//determine what to do

}

template<typename T>
void DriverImaginaryAxis<T>::set_input(input::Setup const& setup) {

}

template<typename T>
bool DriverImaginaryAxis<T>::converged() const {
	return true;
}

} /* namespace eliashberg */
} /* namespace scallop */
