/*	This file EliashbergImaginaryAxis.hpp is part of scallop.
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

#include "scallop/eliashberg/EliashbergImaginaryAxis.h"
#include "scallop/auxillary/BasicFunctions.h"
#include "scallop/auxillary/MixingModule.h"

namespace scallop {
namespace eliashberg {

template<typename T>
EliashbergImaginaryAxis<T>::EliashbergImaginaryAxis(
		MatzubaraEffectiveCouplingMatrix<T> const& diagonalCouplingUpSpin,
		MatzubaraEffectiveCouplingMatrix<T> const& diagonalCouplingDownSpin,
		MatzubaraEffectiveCouplingMatrix<T> const& offDiagonalCoupling) :
				_diagonalCouplingUpSpin(diagonalCouplingUpSpin),
				_diagonalCouplingDownSpin(diagonalCouplingDownSpin),
				_offDiagonalCoupling(offDiagonalCoupling){ };

template<typename T>
void EliashbergImaginaryAxis<T>::solve(
		T temperature,
		EliashbergGapFunction<T> const& deltaEInitalGuess,
		FrequencyCorrectionZ<T> const& ZInitalGuess,
		ASymmetricEnergyCorrection<T> const& AInitialGuess) {

	_deltaE = deltaEInitalGuess;
	auxillary::MixingModule< EliashbergGapFunction<T> > MixingDeltaE(_deltaE);
	_eliashZ = ZInitalGuess;
	auxillary::MixingModule< FrequencyCorrectionZ<T> > MixingZ(_eliashZ);
	_eliashA = AInitialGuess;
	auxillary::MixingModule< ASymmetricEnergyCorrection<T> > MixingA(_eliashA);
	_temperature = temperature;
	T beta = auxillary::BasicFunctions::inverse_temperature(_temperature);

	bool converged;

	size_t numIterations = 0;

	do {

		//Compute square root tmp save
		_gapSquareRootSpinUp.compute(_eliashZ,_deltaE,_eliashA,beta);
		_gapSquareRootSpinDown.compute(_eliashZ,_deltaE,_eliashA,beta);

		//Compute energy integrals
		_spinUpM.compute(_gapSquareRootSpinUp,_gapSquareRootSpinDown,_eliashZ,_eliashA,beta);
		_spinDownM.compute(_gapSquareRootSpinUp,_gapSquareRootSpinDown,_eliashZ,_eliashA,beta);
		_eliashN.compute(_gapSquareRootSpinUp,_gapSquareRootSpinDown,_eliashZ,_eliashA);

		//Main Eliashberg, be careful to use the variables from the last iteration consistently
		_deltaE.compute(_offDiagonalCoupling,_eliashN,MixingZ.get_current_iteration(),beta);
		_eliashA.compute(_diagonalCouplingUpSpin,_diagonalCouplingDownSpin,
				_spinUpM,_spinDownM,MixingZ.get_current_iteration(),beta);
		_eliashZ.compute(_diagonalCouplingUpSpin,_diagonalCouplingDownSpin,
				_spinUpM,_spinDownM,MixingZ.get_current_iteration(),beta);

		//check convergence
		converged = true;
		bool switchToBroydboolen = true;

		//mix in old solutions
		MixingDeltaE.mixing(_deltaE);
		MixingZ.mixing(_eliashZ);
		MixingA.mixing(_eliashA);

		++numIterations;
	} while (not converged);
}

template<typename T>
void EliashbergImaginaryAxis<T>::report_statistics(size_t numIterations, bool converged) {

}

template<typename T>
EliashbergGapFunction<T> const& EliashbergImaginaryAxis<T>::get_eliashberg_gap() const {
	return _deltaE;
}

template<typename T>
FrequencyCorrectionZ<T> const& EliashbergImaginaryAxis<T>::get_eliashberg_Z() const {
	return _eliashZ;
}

template<typename T>
ASymmetricEnergyCorrection<T> const& EliashbergImaginaryAxis<T>::get_eliashberg_A() const {
	return _eliashA;
}

} /* namespace eliashberg */
} /* namespace scallop */
