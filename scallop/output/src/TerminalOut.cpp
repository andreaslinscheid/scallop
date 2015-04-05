/*	This file TerminalOut.cpp is part of scallop.
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
 *  Created on: Dec 5, 2014
 *      Author: Andreas Linscheid
 */

#include "scallop/output/TerminalOut.h"
#include "scallop/auxillary/BasicFunctions.h"
#include "scallop/auxillary/globals.h"
#include <iostream>

namespace scallop {
namespace output {

TerminalOut::TerminalOut() :
		_printToStdErr(false),
		_verbosityLvl(0),
		_buffer(),
		_sstrBuff(_buffer)  { };

TerminalOut::TerminalOut(int verbLvl) :
		_printToStdErr(false),
		_verbosityLvl(verbLvl),
		_buffer(),
		_sstrBuff(_buffer) { };

void TerminalOut::print() const {
	if ( _sstrBuff.str().empty() )
		return;
	if ( auxillary::globals::verbosityLvl < _verbosityLvl )
		return;
	if ( _printToStdErr ){
		std::cerr << _sstrBuff;
	} else {
		std::cout << _sstrBuff;
	}
}

TerminalOut::~TerminalOut() {
	this->print();
}
} /* namespace output */
} /* namespace scallop */
