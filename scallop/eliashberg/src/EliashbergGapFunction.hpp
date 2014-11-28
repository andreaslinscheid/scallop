/*	This file EliashbergGapFunction.hpp is part of scallop.
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

#include "scallop/eliashberg/EliashbergGapFunction.h"

namespace scallop {
namespace eliashberg {

template<typename T>
template<class FreqencyRenormalization,
		 class NambuOffDiagonalEnergyIntegral>
void EliashbergGapFunction<T>::compute(
			T inverseTemperature,
			EliashbergGapFunction<T> const& previousGap,
			MatzubaraEffectiveCouplingMatrix<T> const& couplingNambuDiagonal,
			FreqencyRenormalization const& Z,
			NambuOffDiagonalEnergyIntegral const& N) {

	size_t nMatz = this->get_num_matzubara_pts();
	size_t nSplt = this->get_num_splitting_pts();
	size_t nBand = this->get_num_bands();
	MatzubaraSplittingVector<T> tmp = previousGap;

	for(size_t n=0; n<nMatz ; ++n)
		for(size_t j = 0; j < nSplt ; ++j)
			for(size_t b = 0; b < nBand ; ++b)
				tmp(n,j,b) *= N.read(n,j,b) * Z.read(n,j,b);

	for(size_t n= 0; n < nMatz ; ++n)
		for(size_t j=0 ; j < nSplt ;++j)
			for(size_t b = 0; b < nBand ;++b){

				(*this)(n,j,b) = 0.0;

				for(size_t np=0; np < nMatz ; ++np)
					for(size_t jp=0; jp < nSplt ; ++jp)
						for(size_t bp = 0; bp < nBand ;++bp)
						(*this)(n,j,b) += couplingNambuDiagonal.read(n,j,b,np,jp,bp) * tmp(np,jp,bp);

				(*this)(n,j,b) /= -2.0*inverseTemperature*Z.read(n,j,b);
			}

};

} /* namespace eliashberg */
} /* namespace scallop */
