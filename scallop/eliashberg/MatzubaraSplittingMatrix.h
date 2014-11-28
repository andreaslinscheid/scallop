/*	This file MatzubaraSplittingMatrix.h is part of scallop.
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
 *  Created on: Nov 26, 2014
 *      Author: Andreas Linscheid
 */

#ifndef SCALLOP_ELIASHBERG_MATZUBARASPLITTINGMATRIX_H_
#define SCALLOP_ELIASHBERG_MATZUBARASPLITTINGMATRIX_H_

#include <cstddef>
#include <vector>

namespace scallop {
namespace eliashberg {

template<typename T>
class MatzubaraSplittingMatrix : private std::vector<T> {
public:

	T & operator() (size_t n, size_t j, size_t b, size_t np, size_t jn, size_t bp);

	T read(size_t n, size_t j, size_t b, size_t np, size_t jn, size_t bp) const;

	size_t get_num_matzubara_pts() const;

	size_t get_num_splitting_pts() const;

	size_t get_num_bands() const;
private:

	size_t _matzubaraDim;

	size_t _splittingDim;
};

} /* namespace eliashberg */
} /* namespace scallop */
#include "scallop/eliashberg/src/MatzubaraSplittingMatrix.hpp"
#endif /* SCALLOP_ELIASHBERG_MATZUBARASPLITTINGMATRIX_H_ */
