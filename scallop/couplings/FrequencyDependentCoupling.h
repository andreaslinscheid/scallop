/*	This file FrequencyDependentCoupling.h is part of scallop.
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

#ifndef SCALLOP_COUPLINGS_FREQUENCYDEPENDENTCOUPLING_H_
#define SCALLOP_COUPLINGS_FREQUENCYDEPENDENTCOUPLING_H_

#include "scallop/auxillary/ScallopVector.h"

namespace scallop {
namespace couplings {

template<typename T>
class FrequencyDependentCoupling {

public:
	typedef typename auxillary::ScallopVector<T>::iterator iterator;

	typedef typename auxillary::ScallopVector<T>::const_iterator const_iterator;

	auxillary::ScallopVector<T> const& get_frequency_mesh() const;

	T & operator[] (size_t pos);

	T & operator() (size_t b, size_t j, size_t bp, size_t jp, size_t iomega);

	iterator begin_freq(size_t b, size_t j, size_t bp, size_t jp);
	iterator end_freq(size_t b, size_t j, size_t bp, size_t jp);


private:
	auxillary::ScallopVector<T> _frequencyMesh;
	auxillary::ScallopVector<T> _couplingData;
};

} /* namespace couplings */
} /* namespace scallop */
#include "scallop/couplings/src/FrequencyDependentCoupling.hpp"
#endif /* SCALLOP_COUPLINGS_FREQUENCYDEPENDENTCOUPLING_H_ */
