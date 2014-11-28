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
		MatzubaraSplittingMatrix<double> const& EffectiveCouplingKernelUpSpin,
		MatzubaraSplittingMatrix<double> const& EffectiveCouplingKernelDownSpin,
		MatzubaraSplittingVector<std::complex<double> > const& MUpSpin,
		MatzubaraSplittingVector<std::complex<double> > const& MDownSpin,
		MatzubaraSplittingVector<std::complex<double> > const& Z,
		double inverseTemperature,
		double mixing,
		int digitsToAgreeOnConvergence,
		bool &converged){
	//
	//abbreviation complex identity
	std::complex<double> complex_i = std::complex<double>(0.0,1.0);
	//
	//main computation loop
	for(size_t n=_matzubaraDim;n--;){
		for(size_t j=_splittingDim;j--;){
			//
			//init
			std::complex<double> oldValueForMixing = (*this)(n,j);
			(*this)(n,j) = 0.0;
			//
			//Integration loop
			for(size_t np=_matzubaraDim;np--;){
				for(size_t jp=_splittingDim;jp--;){
					(*this)(n,j) += EffectiveCouplingKernelUpSpin.read(n,j,np,jp)*MUpSpin.read(np,jp);
				}
			}
			for(size_t np=_matzubaraDim;np--;){
				for(size_t jp=_splittingDim;jp--;){
					(*this)(n,j) += EffectiveCouplingKernelDownSpin.read(n,j,np,jp)*MDownSpin.read(np,jp);
				}
			}
			//Scale
			(*this)(n,j) *= complex_i/(4.0*inverseTemperature
					*EliashbergEquations::matzubara_point_to_frequency_fermionic(n-_matzubaraDim/2,inverseTemperature));
			(*this)(n,j) += 1.0;
			//
			//check convergence
			converged = converged and BasicFunctions::agree_on_significant_decimal_digits_with_threshold(
					abs(oldValueForMixing),abs((*this)(n,j)),digitsToAgreeOnConvergence,10e-10);
			//
			//linearly mix with old value
			(*this)(n,j) += (oldValueForMixing - (*this)(n,j))*mixing;
		}
	}
}

} /* namespace eliashberg */
} /* namespace scallop */
