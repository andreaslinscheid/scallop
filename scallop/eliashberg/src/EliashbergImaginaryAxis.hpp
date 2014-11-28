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

namespace scallop {
namespace eliashberg {

template<typename T>
void EliashbergImaginaryAxis<T>::solve(

		T temperature,
		EliashbergGapFunction<T> const& deltaEInitalGuess,
		FrequencyCorrectionZ<T> const& ZInitalGuess) {

	_deltaE = deltaEInitalGuess;
	_Z = ZInitalGuess;
	T beta =

	bool converged = false;

	do {

		//Compute square root tmp save
		_gapSquareRootSpinUp.compute(_Z,_deltaE,_ASymmetricEnergyCorrection,_inverseTemperature,true);
		_gapSquareRootSpinDown.compute(_frequencyCorrection,
				_pairing,_ASymmetricEnergyCorrection,_inverseTemperature,false);
		//
		//Compute energy integrals
		_NambuDiagonalEnergyIntegralUpSpin.compute(_gapSquareRootSpinUp,
				_gapSquareRootSpinDown,_frequencyCorrection,_ASymmetricEnergyCorrection,
				_inverseTemperature,true);
		_NambuDiagonalEnergyIntegralDownSpin.compute(_gapSquareRootSpinUp,
				_gapSquareRootSpinDown,_frequencyCorrection,_ASymmetricEnergyCorrection,
				_inverseTemperature,false);
		_NambuOffDiagonalEnergyIntegral.compute(_gapSquareRootSpinUp,
				_gapSquareRootSpinDown,_frequencyCorrection,_ASymmetricEnergyCorrection);
		//
		//absDiff is the sum of absolut value differences old to new
		bool converged = true;
		bool switchToBroyden = true;
		//
		//Main Eliashberg: compute gap equation
		_pairing.compute(_NambuOffDiagonalEffectiveCoupling,_NambuOffDiagonalEnergyIntegral,
				_frequencyCorrection,_inverseTemperature,_mixingParameter,_convergencyNumberOfDigits,_digitsToSwitchToBroydenMixing,
				_pairingMixing,converged,switchToBroyden);
		//					and Nambu diagonal parts
		_ASymmetricEnergyCorrection.compute(_NambuDiagonalEffectiveCouplingUpSpin,
				_NambuDiagonalEffectiveCouplingDownSpin,_NambuDiagonalEnergyIntegralUpSpin,
				_NambuDiagonalEnergyIntegralDownSpin,_frequencyCorrection,_inverseTemperature,
				_mixingParameter,_convergencyNumberOfDigits,converged);
		//
		_frequencyCorrection.compute(_NambuDiagonalEffectiveCouplingUpSpin,
				_NambuDiagonalEffectiveCouplingDownSpin,_NambuDiagonalEnergyIntegralUpSpin,
				_NambuDiagonalEnergyIntegralDownSpin,_frequencyCorrection,_inverseTemperature,
				_mixingParameter,_convergencyNumberOfDigits,converged);
		//
		//we always consider the loop to be convereged if the gap went below the minimum
		if ( ( fabs(_pairing.get_value_with_max_abs().real()) < fabs(_gapConsideredZero.real()) ) and
				( fabs(_pairing.get_value_with_max_abs().imag()) < fabs(_gapConsideredZero.imag())) )
			converged = true;
		if ( _minNumberOfIterations >= numberOfIterations){
			converged = false;
		}
		//
		//switch to non-linear Broyden Mixing at given convergence
		if ( not (_broydenMixingInUse) and switchToBroyden ){
			std::cout << "\t switched to Broyden Mixing" << "\n";
			restart_mixing(1);
		}
		//
		//report the norm to give the user a hind about the convergence improvement
		std::cout << "\tThe norm of the gap function is " << _pairing.compute_norm_sqrt_of_abs_square() << "\n";
	} while (not converged);
}

} /* namespace eliashberg */
} /* namespace scallop */
