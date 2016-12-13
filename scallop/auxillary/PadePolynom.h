/*	This file PadePolynom.h is part of scallop.
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
 *  Created on: Nov 29, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_AUXILLARY_PADEPOLYNOM_H_
#define SCALLOP_AUXILLARY_PADEPOLYNOM_H_

#include <complex>
#include <vector>

namespace scallop
{
namespace auxillary
{
/**
 * 	Construct an analytic continuation via the Pade method.
 *
 * 	The function must not be initialized for symmetric data w.r.t. the Matsubara frequency.
 * 	Use only one half of the axis instead.
 */
class PadePolynom
{
public:
	typedef std::complex<double> T;

	/**
	 * Initialize the polynom.
	 *
	 * @param coords 	The complex coordinates
	 * @param fvalues	The complex values of the function at \p coords.
	 */
	PadePolynom( std::vector<T> const& coords, std::vector<T> const& fvalues);

	/**
	 * Initialize the polynom.
	 *
	 * @param coords	Point to the complex coordinates
	 * @param fvalues	Point to the complex values of the function at \p coords.
	 * @param dim		The numer of elements in coords and fvalues.
	 */
	PadePolynom( T const * coords, T const * fvalues, size_t dim);

	/**
	 * Evaluate the polynom.
	 *
	 * @param z	The complex number the function should be computed at.
	 * @return	Polynom value at \p z
	 */
	T operator() (T z) const;
private:

	bool isConstant_ = true;

	T constant_ = T(0);

	std::vector<T> coefficients_;
	std::vector<T> initpoints_;

	void compute_polynom( T const * coords, T const * fvalues, size_t dim );
};

} /* namespace auxillary */
} /* namespace scallop */

#endif /* SCALLOP_AUXILLARY_PADEPOLYNOM_H_ */
