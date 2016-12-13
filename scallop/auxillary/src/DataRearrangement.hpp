/*	This file DataRearrangement.hpp is part of scallop.
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
 *  Created on: Dec 1, 2016
 *      Author: A. Linscheid
 */

#include "scallop/auxillary/DataRearrangement.h"
#include <set>

namespace scallop
{
namespace auxillary
{

template<typename T>
void DataRearrangement<T>::redistribute_locally(
		T & data,
		std::vector<size_t> const& gridIndicesMapOldToNew,
		size_t blockSize)
{
	if ( buffer_.size() != blockSize )
		buffer_ = T(blockSize);

	//reshuffel the data according to the index map
	std::set<size_t> mappedIndices;
	std::vector< std::vector<std::pair<size_t,size_t> > > copyLoops;
	for ( size_t i = 0; i < gridIndicesMapOldToNew.size(); ++i )
	{
		if ( mappedIndices.find(i) != mappedIndices.end() )
			continue;
		std::vector<std::pair<size_t,size_t> > loop;
		size_t step = i;
		do
		{
			mappedIndices.insert(step);
			if ( step == gridIndicesMapOldToNew[step] )
				break;
			loop.push_back( std::make_pair(step,gridIndicesMapOldToNew[step]) );
			step = gridIndicesMapOldToNew[step];
		} while ( step != i );
		if ( ! loop.empty() )
			loop.pop_back();//we do not need the last copy that closes the loop

		if ( ! loop.empty() )
			copyLoops.push_back( std::move(loop) );
	}
	for ( auto l : copyLoops )
	{
		//fill buffer with the first target
		std::copy( data.begin()+l.front().first*blockSize,
				data.begin()+(l.front().first+1)*blockSize,
				buffer_.begin());

		for ( auto element : l )
		{
			//reshuffel data. The loop makes sure we never overwrite data
			// that is not already copied inside the array or in the buffer.
			std::copy( data.begin()+element.second*blockSize,
					data.begin()+(element.second+1)*blockSize,
					data.begin()+element.first*blockSize);
		}

		//finally, copy the buffer to the last source
		std::copy(buffer_.begin(),
				buffer_.end(),
				data.begin()+l.back().second*blockSize);
	}
}

} /* namespace auxillary */
} /* namespace scallop */
