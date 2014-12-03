/*	This file GapSquareRoot.hpp is part of scallop.
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
 *  Created on: Nov 27, 2014
 *      Author: Andreas Linscheid
 */

#include "scallop/eliashberg/GapSquareRoot.h"

namespace scallop {
namespace eliashberg {

template<typename T, int spin>
void GapSquareRoot<T,spin>::compute(FrequencyCorrectionZ<T> const& Z,
		EliashbergGapFunction<T> const& deltaE,
		ASymmetricEnergyCorrection<T> const& A,
		T inverseTemperature){

	//abbreviation for the complex identity
	std::complex<double> complex_i = std::complex<double>(0.0,1.0);

	int numMatzPtsDiv2 = static_cast<int>(this->get_num_matzubara_pts())/2;

	//main computation loop
	for (size_t b = 0; b < this->get_num_bands(); ++b )
		for (size_t j = 0; j < this->get_num_splitting_pts(); ++j )
			for (int n = -numMatzPtsDiv2; n < numMatzPtsDiv2; ++n){

				T omegaN = auxillary::BasicFunctions::matzubara_frequency_of_index(n,inverseTemperature);

				(*this)(b,j,n) = std::sqrt	(-( std::pow( Z.read(b,j,n)*deltaE.read(b,j,n) , 2)
											   -std::pow( complex_i * omegaN * Z.read(b,j,n)
														 -spin*_splittingMesh[j], 2)
											  )
										  	);
			}
}
} /* namespace eliashberg */
} /* namespace scallop */
