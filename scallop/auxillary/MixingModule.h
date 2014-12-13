/*	This file MixingModule.h is part of scallop.
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

#ifndef SCALLOP_AUXILLARY_MIXINGMODULE_H_
#define SCALLOP_AUXILLARY_MIXINGMODULE_H_

namespace scallop {
namespace auxillary {

template<class T>
class MixingModule {
public:
	MixingModule(T const& quantityToMix,size_t numberIterationsConsidered = 1);

	void init(T const& quantityToMix,size_t numberIterationsConsidered);
	void restart(size_t numberIterationsConsidered);

	/**
	 *	Apply the modified Broyden mixing.
	 *
	 *	 Ref.	Phys. Rev. B 38, 12807–12813 (1988)
	 *
	 *	On input:
	 *		vectorDimension
	 *
	 *	On output:
	 *		pointerToVectorOld : overwritten with the mixed vector
	 */
	void mixing(double mixingParameter,T & pointerToVectorNewAndOut);

	T const& get_current_iteration() const;

	T const& get_privious_iteration(size_t numItInThePast=1) const;
private:

	std::vector<T> _history;

	///this index mapping is such that _history[_locationConologicalHistory[i]] refers to the ith last added element
	std::vector<size_t> _locationConologicalHistory;

	typename T::value_type _weightsForErrorInverseJacobian;
	std::vector<typename T::value_type> _weightsForErrorIteration;

	std::vector<typename T::value_type> _dotProductBuffer;

	std::vector<typename T::value_type> _beta;
	std::vector<typename T::value_type> _cOfEq15b;
	std::vector<T> _deltaVector;
	std::vector<T> _deltaFunctional;
	bool _isInitialized;
	size_t _numberOfIterationsSoFar;
	size_t _maxNumberOfIterationsConsideredForMixing;
	size_t _vectorDim;

};

} /* namespace auxillary */
} /* namespace scallop */
#include "scallop/auxillary/src/MixingModule.hpp"
#endif /* SCALLOP_AUXILLARY_MIXINGMODULE_H_ */
