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
void MatzubaraEffectiveCouplingMatrix<T>::compute(
		size_t numberMatzubaraPts,
		T temperature,
		couplings::FrequencyDependentCoupling<T> const& couplingData){

	std::vector<T> newContent(numberMatzubaraPts*numberMatzubaraPts,0.0);
	this->initialize(newContent,
			couplingData.get_splitting_mesh(),
			couplingData.get_num_bands(),
			numberMatzubaraPts);

	T const beta = auxillary::BasicFunctions::inverse_temperature(temperature);

	//main computation loop (skipping the upper right )
	for (size_t b=0;b < this->get_num_bands(); ++b)
		for (size_t j=0;j < this->get_num_splitting_pts(); ++j)
			for (size_t n=0;n < this->get_num_matzubara_pts() ; ++n){

				int matzubaraIndex= n - static_cast<int>(this->get_num_matzubara_pts()/2);
				T matzubaraFeqN =
					auxillary::BasicFunctions::matzubara_frequency_of_index(matzubaraIndex,beta);

				for (size_t bp=0;bp < this->get_num_bands(); ++bp)
					for (size_t jp=0;jp < this->get_num_splitting_pts(); ++jp)
						for (size_t np=0;np < n ; ++np){

							int matzubaraIndexP= np - static_cast<int>(this->get_num_matzubara_pts()/2);
							T matzubaraFeqNPrime =
								auxillary::BasicFunctions::matzubara_frequency_of_index(matzubaraIndexP,beta);

							T matzubaraFreqDiffSquare = std::pow(matzubaraFeqN-matzubaraFeqNPrime,2);

							//Integrate phonon coupling
							T frequencyMeshIntervalLength = couplingData.get_frequency_mesh()[1] -
															couplingData.get_frequency_mesh()[0];
							typename couplings::FrequencyDependentCoupling<T>::iterator itCoupl
								= couplingData.begin_freq(b,j,bp,jp);
							for (size_t ifreq = 0; ifreq < couplingData.get_frequency_mesh().size() ; ++ifreq){
								(*this)(b,j,n,bp,jp,np) +=
										2.0*couplingData.get_frequency_mesh()[ifreq]
										   * (*itCoupl)
										/(matzubaraFreqDiffSquare+std::pow(couplingData.get_frequecny_mesh()[ifreq],2))
										*frequencyMeshIntervalLength;
								++itCoupl;
							}
					}
			}

	//employ the symmetry of this object: n <-> n'
	for (size_t b=0;b < this->get_num_bands(); ++b)
		for (size_t j=0;j < this->get_num_splitting_pts(); ++j)
			for (size_t n=this->get_num_matzubara_pts()/2;n < this->get_num_matzubara_pts() ; ++n)
				for (size_t bp=0;bp < this->get_num_bands(); ++bp)
					for (size_t jp=0;jp < this->get_num_splitting_pts(); ++jp)
						for (size_t np=0;np < n ; ++np)
							(*this)(b,j,np,bp,jp,n) = (*this)(b,j,n,bp,jp,np);
}

} /* namespace eliashberg */
} /* namespace scallop */
