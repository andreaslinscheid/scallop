/*	This file Setup.h is part of scallop.
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

#ifndef SCALLOP_INPUT_SETUP_H_
#define SCALLOP_INPUT_SETUP_H_

#include "scallop/input/InputFile.h"
#include <string>
#include <map>

namespace scallop {
namespace input {

/**
 * \brief Perform the initial startup of the code and trigger the InputFile reading.
 *
 * The constructor will fill the internal InputFile with data which can then be accessed by
 * other input classes that parse selected input options from the InputFile.
 * The InputFile can be accessed using \ref get_parsed_input_file()
 */
class Setup {
public:

	/**
	 * \brief Create the initial startup.
	 *
	 * @param argc The number of arguments taken from the input.
	 * @param argv The array of size argc with c-stype char arrays.
	 */
	Setup(int argc, char *argv[]);

	/**
	 * \brief Allows constant access to the internal InputFile.
	 *
	 * Input classes can look up their respective defined input options using the returned const reference.
	 *
	 * @return A const reference to the internal InputFile
	 */
	InputFile const& get_parsed_input_file() const;

	/**
	 * \brief Forward the call to \ref InputFile::get_list_unread_input_parameters()
	 *
	 * @return A list of keys in the input file that have not been referenced.
	 */
	std::vector<std::string> get_list_unread_input_parameters() const;
private:
	InputFile _inputFile;

};

} /* namespace input */
} /* namespace scallop */
#endif /* SCALLOP_INPUT_SETUP_H_ */
