/*	This file FrequencyCorrectionZ.hpp is part of scallop.
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

#include "scallop/eliashberg/FrequencyCorrectionZ.h"

namespace scallop {
namespace eliashberg {

template<typename T>
void FrequencyCorrectionZ<T>::compute(
		MatzubaraEffectiveCouplingMatrix<T> const& EffectiveCouplingKernelUpSpin,
		MatzubaraEffectiveCouplingMatrix<T> const& EffectiveCouplingKernelDownSpin,
		NambuDiagonalEnergyIntegralM<T,auxillary::Constants::upspin> const& MUpSpin,
		NambuDiagonalEnergyIntegralM<T,auxillary::Constants::downspin> const& MDownSpin,
		FrequencyCorrectionZ<T> const& previousZ,
		double inverseTemperature ) {

	std::complex<double> complex_i = std::complex<double>(0.0,1.0);

	for (size_t b = 0 ; b < this->get_num_bands() ; ++b )
		for(size_t j= 0 ; j <  this->get_num_splitting_pts() ; ++j)
			for(size_t n = 0 ; n < this->get_num_matzubara_pts() ; ++n){

				(*this)(b,j,n) = 0.0;

				//Integration loop for up and down spin
				for (size_t bp = 0 ; bp < this->get_num_bands() ; ++bp )
					for(size_t jp= 0 ; jp <  this->get_num_splitting_pts() ; ++jp)
						for(size_t np = 0 ; np < this->get_num_matzubara_pts() ; ++np)
							(*this)(b,j,n) +=
									EffectiveCouplingKernelUpSpin(b,j,n,bp,jp,np)*MUpSpin(bp,jp,np);

				for (size_t bp = 0 ; bp < this->get_num_bands() ; ++bp )
					for(size_t jp= 0 ; jp <  this->get_num_splitting_pts() ; ++jp)
						for(size_t np = 0 ; np < this->get_num_matzubara_pts() ; ++np)
							(*this)(b,j,n) +=
									EffectiveCouplingKernelDownSpin(b,j,n,bp,jp,np)*MDownSpin(bp,jp,np);

				//Scale
				(*this)(b,j,n) *= complex_i/(4.0*inverseTemperature
						*auxillary::BasicFunctions::matzubara_frequency_of_index(n-_matzubaraDim/2,inverseTemperature));
				(*this)(b,j,n) += 1.0;
			}
}

} /* namespace eliashberg */
} /* namespace scallop */
