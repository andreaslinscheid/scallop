/*	This file PadePolynom.cpp is part of scallop.
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

#include "scallop/auxillary/PadePolynom.h"
#include "scallop/error_handling/Error.h"

namespace scallop
{
namespace auxillary
{

PadePolynom::PadePolynom( std::vector<T> const& coords, std::vector<T> const& fvalues )
{
	if ( coords.size() != fvalues.size() )
		error_handling::Error("Incorrect input sizes: coordinate size must match values!",2);

	this->compute_polynom( coords.data(), fvalues.data(), coords.size() );
}

PadePolynom::PadePolynom( T const * coords, T const * fvalues, size_t dim )
{
	this->compute_polynom(coords,fvalues,dim);
}

void PadePolynom::compute_polynom( T const * coords, T const * fvalues, size_t dim )
{
	//The algorithm does not work for a complex constant function. Thus, we specify this case explicitly.
	constant_ = fvalues[0];
	for ( size_t i = 0 ; i < dim ; ++i )
		isConstant_ = isConstant_ && std::abs(fvalues[i]-constant_)<1e-50;

	if ( isConstant_ )
		return;

	initpoints_.assign(coords,coords+dim);

	std::vector<T> tmpRecursionValues( dim*dim, T(0) );
	for ( size_t i = 0; i < dim; i ++)
		tmpRecursionValues[i] = fvalues[i];

	/*
	 * 	Implementation of the recursion formula
	 */
	for (size_t i=1;i<dim;i++)
		for (size_t j=i;j<dim;j++)
		{
			tmpRecursionValues[i*dim+j]=(tmpRecursionValues[(i-1)*dim+i-1]
						-tmpRecursionValues[(i-1)*dim+j])
					/( (initpoints_[j]-initpoints_[i-1])*tmpRecursionValues[(i-1)*dim+j] );

#ifdef DEBUG_BUILD
			if ( tmpRecursionValues[i*dim+j] == T(0) )
			{
				std::cout << "i=" << i << " and j=" << j
				       << " produced a complex zero! ( value one " << tmpRecursionValues[(i-1)*dim+i-1]
				       << " value two" << tmpRecursionValues[(i-1)*dim+j] << ")" <<
				       "\n Furthermore "  << (tmpRecursionValues[(i-1)*dim+i-1]
				                         						-tmpRecursionValues[(i-1)*dim+j]) << std::endl;
				scallop::error_handling::Error(
						std::string("Problem fitting Pade polynom: Recursion value number i=")
					+std::to_string(i)+" and j="+std::to_string(j)+" is a complex zero",2);
			}
#endif
		}

	coefficients_.clear();
	coefficients_.reserve(dim);
	for (size_t i=0;i<dim;i++)
		coefficients_.push_back( tmpRecursionValues[i*dim+i] );
#ifdef DEBUG_BUILD

	for (size_t i=0;i<dim;i++)
		if ( coefficients_[i] != coefficients_[i] )
			scallop::error_handling::Error(std::string("Problem fitting Pade polynom: value number ")
				+std::to_string(i)+" is nan",2);
#endif
}

typename PadePolynom::T PadePolynom::operator() (T z) const
{
	if ( isConstant_ )
		return constant_;

	T previousDenominoator = T(1.0);
	for (size_t i = initpoints_.size()-1 ; i>0 ; i--)
		previousDenominoator = 1.0 + (coefficients_[i]*(z-initpoints_[i-1])/previousDenominoator);

	return coefficients_[0]/previousDenominoator;
}

} /* namespace auxillary */
} /* namespace scallop */
