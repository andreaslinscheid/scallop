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

#include "scallop/auxillary/globals.h"
#include "scallop/parallel/MPIModule.h"
#include <cstddef>
#include <sstream>

namespace scallop
{
namespace output
{

class TerminalOut
{
public:
	TerminalOut();

	TerminalOut( auxillary::globals::VerbosityLvl verbLvl );

	void print();

	template<typename T>
	void insert(T const& msg);

	void print_startup_message();
private:

	const bool printToStdErr_;

	/// Verbosity level when the currently buffered message is printed.
	const auxillary::globals::VerbosityLvl verbosityLvl_;

	std::stringstream sstrBuff_;
};

class MessageChain
{
public:
	MessageChain(TerminalOut & theMessage);

	template<class T>
	void insert(T const& msg) const;

	~MessageChain();
private:
	void operator= (MessageChain rhs) = delete;

	TerminalOut & theMessage_;
};

} /* namespace output */
} /* namespace scallop */
#include "scallop/output/src/TerminalOut.hpp"
#endif /* SCALLOP_OUTPUT_TERMINALOUT_H_ */
