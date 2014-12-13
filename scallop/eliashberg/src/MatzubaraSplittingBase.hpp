/*	This file MatzubaraSplittingBase.hpp is part of scallop.
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
 *  Created on: Nov 28, 2014
 *      Author: Andreas Linscheid
 */

#include "scallop/eliashberg/MatzubaraSplittingBase.h"
#include "scallop/error_handling/Error.h"

namespace scallop {
namespace eliashberg {

template<class derived,typename T>
size_t MatzubaraSplittingBase<derived,T>::get_num_matzubara_pts() const {
	return _numMatzPts;
}

template<class derived,typename T>
size_t MatzubaraSplittingBase<derived,T>::get_num_splitting_pts() const {
	return _numSplitPts;
}

template<class derived,typename T>
size_t MatzubaraSplittingBase<derived,T>::get_num_bands() const {
	return _numBands;
}

template<class derived,typename T>
T & MatzubaraSplittingBase<derived,T>::operator[] (size_t i) {
#ifdef DEBUG_BUILD
	if ( not _isInit )
		error_handling::Error("Access on not initialized object", 2);
	if ( i >= this->size() )
		error_handling::Error("out of range access", 2);
#endif
	return std::vector<T>::operator [](i);
}

template<class derived,typename T>
T MatzubaraSplittingBase<derived,T>::operator[] (size_t i) const{
	return (*this)[i];
}

template<class derived,typename T>
void MatzubaraSplittingBase<derived,T>::initialize(std::vector<T> & newContent,
		std::vector<T> const& splittingPoints,
		size_t NumBands,
		size_t MatzubaraPts) {

	if ( _isInit )
		error_handling::Error("We do not allow to re-initialize MatzubaraSplittingBase", 2);

	_isInit = true;

	_splittingPoints = splittingPoints;

	_numBands = NumBands;
	_numSplitPts = _splittingPoints.size();
	_numMatzPts = MatzubaraPts;

	//slice the vector from this object and replace it by newContent
	std::vector<T> & baseThis = *this;
	std::swap(baseThis,newContent);
}

template<class derived,typename T>
T MatzubaraSplittingBase<derived,T>::get_splitting_pt(size_t index) const {
#ifdef DEBUG_BUILD
	if ( index >= _splittingPoints.size() )
		error_handling::Error("Splitting points out of range access", 2);
#endif
	return _splittingPoints[index];
}

} /* namespace eliashberg */
} /* namespace scallop */
