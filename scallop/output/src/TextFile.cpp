/*	This file TextFile.cpp is part of scallop.
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

#include "scallop/output/TextFile.h"
#include "scallop/error_handling/Error.h"
#include <fstream>
#include <cerrno>
#include <cstring>

namespace scallop {
namespace output {

void TextFile::write(std::string const& fileName,
		std::string const& content)  const {
	std::ofstream file;
	errno = 0;
	file.open(fileName.c_str(),std::ios_base::out);
	if ( not file.good() )
		error_handling::Error(std::string()+"Cannot open file "+fileName+" for reading: "+std::strerror(errno), 1);
	file << content;
}

} /* namespace output */
} /* namespace scallop */
