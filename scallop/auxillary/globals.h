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

#define SCALLOP_MAJOR_VERSION 0
#define SCALLOP_MINOR_VERSION 1

///Set the global verbosity level of command line output.
/// The numbers must be such that high > medium > low
typedef enum
{
	high = 3000,
	medium = 2000,
	low = 1000
} VerbosityLvl;

///Before output, compare if a local verbosity level is <= than this global number and only print if it is.
extern VerbosityLvl vLvl;

} /*namespace globals */
} /*namespace auxillary */
} /*namespace scallop */

#endif /* SCALLOP_AUXILLARY_GLOBALS_H_ */
