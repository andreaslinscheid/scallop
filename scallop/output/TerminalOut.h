/*	This file TerminalOut.h is part of scallop.
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

#ifndef SCALLOP_OUTPUT_TERMINALOUT_H_
#define SCALLOP_OUTPUT_TERMINALOUT_H_

#include <cstddef>
#include <sstream>

namespace scallop {
namespace output {

class TerminalOut {
public:

	enum verbosityLvl {
		high = 1000,
		medium = 1100,
		low = 1110
	};

	TerminalOut();

	~TerminalOut();

	void print() const;

	TerminalOut(int verbLvl);

	/**
	 * \brief Allow the user to insert any type of message.
	 *
	 * It can be used with any type that can be stringified using the operator<< of a stringstream
	 *
	 * @param message The message the internal buffer is set to.
	 */
	template<typename T>
	TerminalOut& operator<< (T const& message);

private:

	const bool _printToStdErr;

	/// Verbosity level when the currently buffered message is printed.
	/// 0 means low the message is always printed, while 1 will only be printed if
	/// the global verbosity level is higher or equal than medium.
	const int _verbosityLvl;

	std::string _buffer;

	std::stringstream _sstrBuff;
};

} /* namespace output */
} /* namespace scallop */
#include "scallop/output/src/TerminalOut.hpp"
#endif /* SCALLOP_OUTPUT_TERMINALOUT_H_ */
