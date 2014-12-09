/*	This file MixingModule.hpp is part of scallop.
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
 *  Created on: Nov 30, 2014
 *      Author: Andreas Linscheid
 */
#include "scallop/auxillary/MixingModule.h"
#include "scallop/error_handling/Error.h"

namespace scallop {
namespace auxillary {

template<class T>
MixingModule<T>::MixingModule(T const& quantityToMix,size_t numberIterationsConsidered)
	:	_history(numberIterationsConsidered,quantityToMix), _refCurrentIteration(&_history[0]) { };


template<class T>
void MixingModule<T>::init(T const& quantityToMix,size_t numberIterationsConsidered) {
	if ( numberIterationsConsidered > 0)
		error_handling::Error("Mixing must use 1 or more previous steps",2);
	_history = std::vector<T>(numberIterationsConsidered,quantityToMix);
	_refStepsHistory.clear();
	_refStepsHistory.reserve(numberIterationsConsidered);
	for ( auto &&refH : _history)
		_refStepsHistory.push_back(refH);
}

template<class T>
void MixingModule<T>::mixing(double mixingParameter,T & pointerToVectorNewAndOut) {
	for ( auto it = _history.begin()+1 ; it != _history.end(); ++it)
		*it = *(it - 1);

	pointerToVectorNewAndOut
}

} /* namespace auxillary */
} /* namespace scallop */
