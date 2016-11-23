/*	This file InteractionMatrix.hpp is part of scallop.
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
 *  Created on: Nov 22, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/InteractionMatrix.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
void InteractionMatrix<T>::init_file( std::string const& filename )
{
	std::ifstream file( filename.c_str() );
	if ( ! file.good() )
		error_handling::Error( std::string("Problem opening file ")+filename,5);

	struct matel
	{
		matel() {};

		matel(size_t j,size_t jp, size_t l1, size_t l2, size_t l3, size_t l4, bT vRe, bT vIm )
		{
			size_t a[6] = {j,jp,l1,l2,l3,l4};
			std::copy(a + 0, a + 6, ind);
			m = T(vRe,vIm);
		};

	    size_t ind[6];
	    T m = T(0);

	    bool operator< ( matel const& other ) const
	    {
	    	for (size_t i = 0 ; i < 5; ++i)
	    		if (ind[i] != other.ind[i] )
	    			return ind[i]<other.ind[i];
	    	return false;
	    };
	};

	std::string line;
	std::getline(file, line);
    std::istringstream ish(line);
    size_t nOrb;
    ish >> nOrb;

	std::set<matel> elements;
	while (std::getline(file, line))
	{
	    std::istringstream iss(line);
	    matel m;
    	for (size_t i = 0 ; i < 6; ++i)
    		if (!(iss >> m.ind[i]))
    			scallop::error_handling::Error( std::string("Expected a 8 column file but found line ")+line,5);
    	bT re,im;
		if (!(iss >> re ))
			scallop::error_handling::Error( std::string("Expected a 8 column file but found line ")+line,5);
		if (!(iss >> im ))
			scallop::error_handling::Error( std::string("Expected a 8 column file but found line ")+line,5);
		m.m = T(re,im);

		auto it = elements.insert( std::move(m) );
		if ( ! it.second )
			scallop::error_handling::Error( std::string("Interaction matrix element defined twice in file ")+filename,5);
	}

	this->initialize_layout_4pt_scalar_obj( nOrb, 4 );

	data_.assign( 16*std::pow(nOrb,4), T(0) );
	for (size_t j = 0 ; j < 4; j++)
		for (size_t jp = 0 ; jp < 4; jp++)
			for (size_t l1 = 0 ; l1 < nOrb; l1++)
				for (size_t l2 = 0 ; l2 < nOrb; l2++)
					for (size_t l3 = 0 ; l3 < nOrb; l3++)
						for (size_t l4 = 0 ; l4 < nOrb; l4++)
						{
							auto it = elements.find( matel(j,jp,l1,l2,l3,l4,bT(0),bT(0)) );
							(*this)(j,jp,l1,l2,l3,l4) = ( it != elements.end() ? it->m : T(0) );
						}
}

template<typename T>
T & InteractionMatrix<T>::operator() (size_t j, size_t jp, size_t l1, size_t l2, size_t l3, size_t l4)
{
	return *(this->write_ptr() + this->memory_layout_4pt_scalar_obj(j,jp,l1,l2,l3,l4));
}

template<typename T>
T * InteractionMatrix<T>::write_ptr()
{
	return data_.data();
}

template<typename T>
T const * InteractionMatrix<T>::read_ptr() const
{
	return data_.data();
}

} /* namespace gw_flex */
} /* namespace scallop */
