/*	This file Setup.cpp is part of scallop.
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
 *  Created on: Nov 14, 2014
 *      Author: Andreas Linscheid
 */

#include "scallop/input/Setup.h"
#include "scallop/error_handling/Error.h"
#include "scallop/output/TerminalOut.h"
#include <fstream>
#include <vector>
#include <iostream>

namespace scallop {
namespace input {

Setup::Setup(int argc, char *argv[])
	:  _inputFile(), config_()
{
	//If we are in debugging mode, generate a stack trace in the case of a segmentation fault
	error_handling::Error::generate_stacktrace_on_segfault();

	//Print header
	output::TerminalOut msg;
	msg.print_startup_message();

	//Read the input file
	if ( argc < 1 ) {
		error_handling::Error( "\tNo options passed, exiting ... !", 1 );
	}else if ( argc == 2 ) {
		//fill the internal key/value map with content.
		_inputFile.read_input_file(argv[1]);
	} else {
		error_handling::Error( "\tCannot figure out how to get input, exiting ... !", 1 );
	}

	//Fetch the general configuration variables from the file
	config_.parse_variables(_inputFile);

	//Flag a warning for all keys that have not been used, as this indicates unintended input
	std::vector<std::string> const unusedKeys = this->get_list_unread_input_parameters();
	for ( auto & key : unusedKeys )
		error_handling::Warning("Key '"+ key +"' defined in the input file has not been\n"
				" matched to a variable used by the code!") ;

	if ( config_.get_vbl().compare("l") == 0 )
	{
		auxillary::globals::vLvl = auxillary::globals::VerbosityLvl::low;
	}
	else if (config_.get_vbl().compare("m") == 0)
	{
		auxillary::globals::vLvl = auxillary::globals::VerbosityLvl::medium;
	}
	else if (config_.get_vbl().compare("h") == 0)
	{
		auxillary::globals::vLvl = auxillary::globals::VerbosityLvl::high;
	}
	else
	{
		error_handling::Error("Unkown vbl value!");
	}

	//In the following, we perform input checks if the input is logically sound
	if ( (not config_.get_r_omega().empty() ) && (config_.get_r_omega().size() != 4) )
		error_handling::Error("r_omega needs 4 values: min, max, #pts, and complex eta!");

	if ( config_.get_f_kpath().empty() && ( (not config_.get_f_spec().empty()) ) )
		error_handling::Error("Output requires the f_kpath to be set to a valid k path file.");
}


InputFile const& Setup::get_parsed_input_file() const
{
	return _inputFile;
}

std::vector<std::string> Setup::get_list_unread_input_parameters() const
{
	return _inputFile.get_list_unread_input_parameters();
}

Configuration const& Setup::get_config() const
{
	return config_;
}

} /* namespace input */
} /* namespace scallop */
