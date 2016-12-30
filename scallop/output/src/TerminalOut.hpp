/*	This file TerminalOut.hpp is part of scallop.
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

namespace scallop
{
namespace output
{

template<class T>
MessageChain operator<<(TerminalOut &msg, T const data)
{
	msg.insert(data);
	return MessageChain(msg);
};

template<class T>
MessageChain const& operator<<(MessageChain const& msgChain, T const data)
{
	msgChain.insert(data);
	return msgChain;
};

template<typename T>
void TerminalOut::insert(T const& msg)
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	if ( ! mpi.ioproc() )
		return;

	if ( ! ( verbosityLvl_<= auxillary::globals::vLvl  ) )
		return;

	sstrBuff_ << msg;
}

template<class T>
void MessageChain::insert(T const& data) const
{
	theMessage_.insert(data);
}
} /* namespace output */
} /* namespace scallop */
