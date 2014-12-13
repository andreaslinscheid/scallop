/*	This file MatzubaraEffectiveCouplingMatrix.hpp is part of scallop.
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

#include "scallop/eliashberg/MatzubaraEffectiveCouplingMatrix.h"
#include "scallop/auxillary/BasicFunctions.h"
#include <cmath>

namespace scallop {
namespace eliashberg {

template<typename T>
T & MatzubaraEffectiveCouplingMatrix<T>::operator() (size_t b, size_t j, size_t bp, size_t jp, int nMinusNp) {
	return MatzubaraSplittingVector<T>::operator()(b,j,bp,jp,nMinusNp);
}

template<typename T>
T MatzubaraEffectiveCouplingMatrix<T>::operator() (size_t b, size_t j, size_t bp, size_t jp, int nMinusNp) const{
	return MatzubaraEffectiveCouplingMatrix<T>::operator()(b,j,bp,jp,nMinusNp);
}

template<typename T>
void MatzubaraEffectiveCouplingMatrix<T>::compute(
		size_t numberMatzubaraPts,
		T temperature,
		couplings::FrequencyDependentCoupling<T> const& couplingData){

	//Initialize itself to zero
	std::vector<T> newContent(numberMatzubaraPts*2*
			couplingData.get_num_splitting_pts()*
			couplingData.get_num_bands(),0.0);
	this->initialize(newContent,
			couplingData.get_splitting_mesh(),
			couplingData.get_num_bands(),
			numberMatzubaraPts*2);

	T const beta = auxillary::BasicFunctions::inverse_temperature(temperature);

	std::vector<T> omegaNMinusOmegaNprimeSqr(numberMatzubaraPts*2);
	for ( size_t n = 0 ; n < omegaNMinusOmegaNprimeSqr.size(); ++n){
		T diffMatzubaraFeqNMinusNPrime = std::pow(M_PI*2*n/beta,2);
		omegaNMinusOmegaNprimeSqr[n] = diffMatzubaraFeqNMinusNPrime;
	}

	//main computation loop (skipping the upper right )
	for (size_t b=0;b < this->get_num_bands(); ++b)
		for (size_t j=0;j < this->get_num_splitting_pts(); ++j)
			for (size_t bp=0;bp < this->get_num_bands(); ++bp)
				for (size_t jp=0;jp < this->get_num_splitting_pts(); ++jp)
					for (size_t nMnp=0;nMnp < omegaNMinusOmegaNprimeSqr.size() ; ++nMnp){

						//Integrate phonon coupling
						T frequencyMeshIntervalLength = couplingData.get_frequency_mesh()[1] -
														couplingData.get_frequency_mesh()[0];
						typename couplings::FrequencyDependentCoupling<T>::iterator itCoupl
							= couplingData.begin_freq(b,j,bp,jp);
						for (size_t ifreq = 0; ifreq < couplingData.get_frequency_mesh().size() ; ++ifreq){
							(*this)(b,j,n,jp,nMnp-this->get_num_matzubara_pts()/2) +=
									2.0*couplingData.get_frequency_mesh()[ifreq]
									   * (*itCoupl)
									/(omegaNMinusOmegaNprimeSqr[nMnp]+std::pow(couplingData.get_frequecny_mesh()[ifreq],2))
									*frequencyMeshIntervalLength;
							++itCoupl;
						}
					}
}

} /* namespace eliashberg */
} /* namespace scallop */
