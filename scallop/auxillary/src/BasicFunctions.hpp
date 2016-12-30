/*	This file BasicFunctions.hpp is part of scallop.
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

#include "scallop/auxillary/BasicFunctions.h"
#include <cmath>
#include <assert.h>
#include <iostream>

namespace scallop {
namespace auxillary {

template<typename T>
T BasicFunctions::fermi_funcs(T beta, T energy)
{
	return T(1.0)/(std::exp( beta*energy ) + T(1.0) );
}

template<typename T>
T BasicFunctions::inverse_temperature(T temperature){
	return 1.0/(Constants<T>::kBoltzmannMEV*temperature);
}

template<typename T>
std::complex<T>
BasicFunctions::matzubara_frequency_of_index(size_t i, size_t nM, T beta) {
	int frequencyIndex =
			(i < nM/2 ? static_cast<int>(i) : static_cast<int>(i)-static_cast<int>(nM) );
	return std::complex<T>(0,M_PI*(2*frequencyIndex+1)/beta);
}

template<typename T>
typename TemplateTypedefs< std::complex<T> >::scallop_vector
BasicFunctions::matzubara_frequency_array( size_t nM, T beta)
{
	typename TemplateTypedefs< std::complex<T> >::scallop_vector result( nM );
	for ( size_t n = 0 ; n < nM; ++n)
		result[n] = BasicFunctions::matzubara_frequency_of_index(n,nM,beta);
	return result;
}

template<class F, typename T>
void BasicFunctions::secant_method(
		T initial, std::pair<T,T> range,
		F const& function,
		size_t nIterMax,
		std::pair<bool,T> & result) const
{
	size_t nIter = 0;

	T xnm2 = initial;
	T fxnm2 = function(xnm2);
	T xnm1 = initial+std::abs(range.first-range.second)*T(0.01);
	T fxnm1 = function(xnm1);

	result.first = false;
	while ( true )
	{
		result.second = xnm1 - fxnm1 * ( xnm1 - xnm2 )/( fxnm1 - fxnm2 );

		//ensure that we are not leaving the range of the band structure which would be pointless
		if ( result.second <= range.first)
		{
			std::cout << result.second << '\t' << xnm1 << '\t' <<range.first << std::endl;
			result.second = (xnm1+range.first)/T(2);
		}
		if ( result.second >= range.second)
		{
			std::cout << result.second << '\t' << xnm1 << '\t' <<range.second << std::endl;
			result.second = (xnm1+range.second)/T(2);
		}

		T fxn = function(result.second);

		if ( function.check_conv(fxnm1,fxn) )
		{
			result.first = true;
			break;
		}
		xnm2 = xnm1;
		xnm1 = result.second;
		fxnm2 = fxnm1;
		fxnm1 = fxn;

		if ( nIter == nIterMax )
			break;
		nIter++;
	}
}

template<class L, class C, typename T>
void BasicFunctions::secant_method(
		T initial, std::pair<T,T> range,
		L const& function,
		C const& convergenceCheck,
		size_t nIterMax,
		std::pair<bool,T> & result) const
{
	struct Wrapper
	{
		Wrapper(L l, C c) : l_(l), c_(c) { };

		T operator() (T x) const
		{
			return l_(x);
		};

		bool check_conv(T fxn, T fxn_plus_1) const
		{
			return c_(fxn,fxn_plus_1);
		};

		L l_;
		C c_;
	};

	Wrapper funcWrap(function,convergenceCheck);

	this->secant_method(initial,range,funcWrap,nIterMax,result);
}

template<typename T>
typename auxillary::TemplateTypedefs<T>::NambuSpinPauliMatrix
BasicFunctions::get_v(size_t i)
{
	assert(i<4);
	if ( i == 0 )
	{
		auto v1 = [] (size_t a1, size_t s1, size_t a2, size_t s2){
			if ( (a1 == a2) && (s1 == s2 ) )
				return T( a1==0? 1.0 : -1.0 );
			return T(0);
		};
		return v1;
	}
	else if ( i == 1 )
	{
		auto v2 = [] (size_t a1, size_t s1, size_t a2, size_t s2){
			if ( (a1 == a2) && (s1 != s2 ) )
				return T( a1==0? 1.0 : -1.0 );
			return T(0);
		};
		return v2;
	}
	else if ( i == 2 )
	{
		auto v3 = [] (size_t a1, size_t s1, size_t a2, size_t s2){
			if ( (a1 == a2) && (s1 != s2 ) )
				return T(0, s1==0? -1.0 : 1.0 );
			return T(0);
		};
		return v3;
	}

	auto v4 = [] (size_t a1, size_t s1, size_t a2, size_t s2){
		if ( (a1 == a2) && (s1 == s2 ) )
			return T( s1+a1==1? -1.0 : 1.0 );
		return T(0);
	};
	return v4;
}

} /* namespace auxillary */
} /* namespace scallop */
