/*	This file ObservableStatistics.h is part of scallop.
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
 *  Created on: Dec 20, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_OUTPUT_OBSERVABLESTATISTICS_H_
#define SCALLOP_OUTPUT_OBSERVABLESTATISTICS_H_

namespace scallop
{
namespace output
{

template<typename T>
class ObservableStatistics
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;

	void print_statistics( output::TerminalOut & msg,
			gw_flex::SelfEnergy<T> const& Sigma,
			gw_flex::KohnShamBandStructure<T> const& KS,
			gw_flex::GreensFunctionOrbital<T> const& G,
			bT beta) const;
};

} /* namespace output */
} /* namespace scallop */

#include "scallop/output/src/ObservableStatistics.hpp"
#endif /* SCALLOP_OUTPUT_OBSERVABLESTATISTICS_H_ */
