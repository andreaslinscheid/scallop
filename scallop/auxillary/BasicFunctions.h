/*	This file BasicFunctions.h is part of scallop.
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

#ifndef SCALLOP_AUXILLARY_BASICFUNCTIONS_H_
#define SCALLOP_AUXILLARY_BASICFUNCTIONS_H_

#include "scallop/auxillary/TemplateTypedefs.h"
#include <string>
#include <complex>

namespace scallop
{
namespace auxillary
{

struct BasicFunctions {

	template<typename T>
	static T fermi_funcs(T beta, T energy);

	template<typename T>
	static T inverse_temperature(T temperature);

	template<typename T>
	static std::complex<T> matzubara_frequency_of_index(size_t i, size_t nM, T beta);

	template<typename T>
	static typename TemplateTypedefs< std::complex<T> >::scallop_vector
	matzubara_frequency_array( size_t nM, T beta);

	std::string get_scallop_path() const;

	/**
	 * Solve x* = F(x*) via the secand method for x in a given range.
	 *
	 * @param initial		The starting guess for x.
	 * @param range			Boundary of the range to be searched.
	 * @param function		Must implement the method ' T operator() (T xn ) const' to evaluate
	 * 						F and 'bool check_conv(T F_xn, T F_xn_plus_1) const'.
	 * @param nIterMax		Break after this number of iterations.
	 * @param result		Pair where the first parameter indicates success and the second the best approximate to x*.
	 */
	template<class F, typename T>
	void secant_method(
			T initial, std::pair<T,T> range,
			F const& function,
			size_t nIterMax,
			std::pair<bool,T> & result) const;

	/**
	 * Solve x* = L(x*) via the secand method for x in a given range.
	 *
	 * @param initial		The starting guess for x.
	 * @param range			Boundary of the range to be searched.
	 * @param function		Must implement the method ' T operator() (T xn ) const' to evaluate
 	 * @param convergenceCheck Must be functor of signature 'bool convergenceCheck(T F_xn, T F_xn_plus_1) const'.
	 * @param nIterMax		Break after this number of iterations.
	 * @param result		Pair where the first parameter indicates success and the second the best approximate to x*.
	 */
	template<class L, class C, typename T>
	void secant_method(
			T initial, std::pair<T,T> range,
			L const& function,
			C const& convergenceCheck,
			size_t nIterMax,
			std::pair<bool,T> & result) const;

	/**
	 * Return a function point to a Nambu/Spin outer product of Pauli matrices.
	 *
	 * @param i		Index of the particular combination v(i+1=1,...,4)
	 * @return		Function pointer to v(i+1=1,...,4)
	 */
	template<typename T>
	static typename auxillary::TemplateTypedefs<T>::NambuSpinPauliMatrix
	get_v(size_t i);
};

template<typename T>
struct Constants {

	static const int upspin = 1;

	static const int downspin = -1;

	static constexpr T kBoltzmannMHartree = 0.003166811429;

	static constexpr T kBoltzmannMEV = 8.6173324*1e-2;

	static constexpr T THZToMHartree = 0.151983;

	static constexpr T RYToMHartree = 500;

	static constexpr T cmToTheMinus1ToMHartree = 0.0000045563352812*1000;


};

} /* namespace auxillary */
} /* namespace scallop */
#include "scallop/auxillary/src/BasicFunctions.hpp"
#endif /* SCALLOP_AUXILLARY_BASICFUNCTIONS_H_ */
