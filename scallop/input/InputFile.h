/*	This file InputFile.h is part of scallop.
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
 *  Created on: Nov 16, 2014
 *      Author: Andreas Linscheid
 */

#ifndef SCALLOP_INPUT_INPUTFILE_H_
#define SCALLOP_INPUT_INPUTFILE_H_

#include <string>
#include <map>

namespace scallop {
namespace input {

class InputFile {
public:
	/**
	 * Inspect the stored key/value pairs.
	 *
	 * @param key String with the key name.
	 * @return	The value string, or an empty string if there is no such key.
	 */
	std::string get_input_config_value( std::string const& key ) const;

	void read_input_file(std::string const& fileName, std::string &infileContent) const;

	void parse_input(std::string const&input);

	void trim_string(std::string & str, std::string const& whitespace = " \t") const;
private:

	std::map<std::string,std::string> _inputFileKeyValue;
};

} /* namespace input */
} /* namespace scallop */
#endif /* SCALLOP_INPUT_INPUTFILE_H_ */
