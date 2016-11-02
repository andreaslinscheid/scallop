/*	This file TypeMapComplex.h is part of scallop.
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
 *  Created on: Nov 1, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_AUXILLARY_TYPEMAPCOMPLEX_H_
#define SCALLOP_AUXILLARY_TYPEMAPCOMPLEX_H_

#include <complex>

namespace scallop
{
namespace auxillary
{

template<class T>
struct TypeMapComplex
{
	typedef T type;
};

template<class T>
struct TypeMapComplex<std::complex<T> >
{
	typedef T type;
};

}; /* namespace scallop */
}; /* namespace auxillary */

#endif /* SCALLOP_AUXILLARY_TYPEMAPCOMPLEX_H_ */
