/*	This file NambuOffDiagonalEnergyIntegralN.hpp is part of scallop.
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

#include "scallop/eliashberg/NambuOffDiagonalEnergyIntegralN.h"
#include <complex>

namespace scallop {
namespace eliashberg {

template<typename T>
template<class gapSquareRoot, class FrequencyRenorm, class AsymetricEnergyRenorm>
void NambuOffDiagonalEnergyIntegralN<T>::compute (
		gapSquareRoot const& gapSqrtUp,
		gapSquareRoot const& gapSqrtDown,
		FrequencyRenorm const& Z,
		AsymetricEnergyRenorm const& A){

	std::complex<T> iPi(0.0,M_PI);

	for (size_t n=0;n< this->ge;++n)
		for (size_t j=0;j<_splittingDim;++j)
			for (size_t b=0;b<_matzubaraDim;++b){
			//init
			std::complex<double> tmp =0;
			//
			//implement the theta functions
			if ( std::imag(-A.read(n,j)*Z.read(n,j) - gapSqrtUp.read(n,j)) > 0 ){
				tmp += iPi/gapSqrtUp.read(n,j);
			}
			if ( std::imag(-A.read(n,j)*Z.read(n,j) + gapSqrtUp.read(n,j)) > 0 ){
				tmp += -iPi/gapSqrtUp.read(n,j);
			}
			if ( std::imag(+A.read(n,j)*Z.read(n,j) - gapSqrtDown.read(n,j)) > 0 ){
				tmp += iPi/gapSqrtDown.read(n,j);
			}
			if ( std::imag(+A.read(n,j)*Z.read(n,j) + gapSqrtDown.read(n,j)) > 0 ){
				tmp += -iPi/gapSqrtDown.read(n,j);
			}
			(*this)(n,j) = tmp;
		}
}

} /* namespace eliashberg */
} /* namespace scallop */
