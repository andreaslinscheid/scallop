/*	This file TextFile.h is part of scallop.
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
 *  Created on: Nov 18, 2014
 *      Author: Andreas Linscheid
 */

#ifndef SCALLOP_OUTPUT_TEXTFILE_H_
#define SCALLOP_OUTPUT_TEXTFILE_H_

#include <string>

namespace scallop {
namespace output {

/**
 * 	\brief Manage the textFile output.
 *
 * 	Text files are all ASCII.
 */
class TextFile {
public:

	/**
	 * \brief Create or overwrite a file and fill it with content.
	 *
	 * @param fileName The name of the file.
	 * @param content The new content of the file on return.
	 */
	void write(std::string const& fileName, std::string const& content) const;
};

} /* namespace output */
} /* namespace scallop */
#endif /* SCALLOP_OUTPUT_TEXTFILE_H_ */
