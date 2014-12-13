/*	This file MatzubaraSplittingBase.h is part of scallop.
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

#ifndef SCALLOP_ELIASHBERG_MATZUBARASPLITTINGBASE_H_
#define SCALLOP_ELIASHBERG_MATZUBARASPLITTINGBASE_H_

#include <vector>
#include <cstddef>

namespace scallop {
namespace eliashberg {

template<class derived,typename T>
class MatzubaraSplittingBase : private std::vector<T> {

public:

	typedef typename std::vector<T>::iterator iterator;
	typedef typename std::vector<T>::const_iterator const_iterator;

	using std::vector<T>::size;
	using std::vector<T>::begin;
	using std::vector<T>::end;

	const_iterator get_const_iterator(size_t pos) const;

	iterator get_iterator(size_t pos) const;

	T & operator[] (size_t i);

	T operator[] (size_t i) const;

	size_t get_num_matzubara_pts() const;

	size_t get_num_splitting_pts() const;

	size_t get_num_bands() const;

	void initialize(std::vector<T> & newContent,
			std::vector<T> const& splittingPoints,
			size_t NumBands,
			size_t MatzubaraPts);

	T get_splitting_pt(size_t index) const;

private:

	bool _isInit = false;

	size_t _numMatzPts = 0;

	size_t _numSplitPts = 0;

	size_t _numBands = 0;

	std::vector<T> _splittingPoints;
};

} /* namespace eliashberg */
} /* namespace scallop */
#include "scallop/eliashberg/src/MatzubaraSplittingBase.hpp"
#endif /* SCALLOP_ELIASHBERG_MATZUBARASPLITTINGBASE_H_ */
