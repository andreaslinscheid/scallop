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
#include <iostream>

namespace scallop
{
namespace output
{

TerminalOut::TerminalOut() :
		printToStdErr_(false),
		verbosityLvl_(auxillary::globals::VerbosityLvl::high),
		sstrBuff_( std::string() )
{

};

TerminalOut::TerminalOut(auxillary::globals::VerbosityLvl verbLvl) :
		printToStdErr_(false),
		verbosityLvl_(verbLvl),
		sstrBuff_( std::string() )
{

};

void TerminalOut::print()
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();

	if ( mpi.get_mpi_me() != 0 )
		return;

	if ( ! printToStdErr_ )
	{
		std::cout << sstrBuff_.str() << std::endl;
	}
	else
	{
		std::cerr << sstrBuff_.str() << std::endl;
	}
	sstrBuff_.str( std::string() );
	sstrBuff_.clear();
}

MessageChain::~MessageChain()
{
	theMessage_.print();
}


} /* namespace output */
} /* namespace scallop */
