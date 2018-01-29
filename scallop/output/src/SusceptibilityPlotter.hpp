/*	This file SusceptibilityPlotter.hpp is part of scallop.
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
 *  Created on: Feb 7, 2017
 *      Author: A. Linscheid
 */

#include "scallop/output/SusceptibilityPlotter.h"

namespace scallop
{
namespace output
{

template<class D>
void SusceptibilityPlotter::append_data_to_buffer(D const& data)
{
	assert( buffer_.size() >= pos_+data.size() );
	std::copy( data.begin(), data.end(),buffer_.begin()+pos_);
	pos_ += data.size();
}

} /* namespace output */
} /* namespace scallop */
