/*	This file globals.h is part of scallop.
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
 *  Created on: Apr 5, 2015
 *      Author: alinsch
 */

#ifndef SCALLOP_AUXILLARY_GLOBALS_H_
#define SCALLOP_AUXILLARY_GLOBALS_H_

#include <vector>

namespace scallop
{
namespace auxillary
{
namespace globals
{

///Set the global verbosity level of command line output.
typedef enum
{
	high = 1000,
	medium = 1100,
	low = 1110
} VerbosityLvl;

extern VerbosityLvl vLvl;

} /*namespace globals */
} /*namespace auxillary */
} /*namespace scallop */

#endif /* SCALLOP_AUXILLARY_GLOBALS_H_ */
