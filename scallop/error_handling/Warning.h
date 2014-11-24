/*	This file Warning.h is part of scallop.
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
 *  Created on: Nov 24, 2014
 *      Author: Andreas Linscheid
 */

#ifndef SCALLOP_ERROR_HANDLING_WARNING_H_
#define SCALLOP_ERROR_HANDLING_WARNING_H_

namespace scallop {
namespace error_handling {

/**
 * 	\brief Collect Warnings and print them to std::cerr.
 *
 * 	The class collects warnings into an internal buffer that is send
 * 	to std::cerr when it goes out of scope.
 */
class Warning {
public:

	/**
	 * 	\brief Initialize the warning and the empty internal buffer.
	 */
	Warning();

	/**
	 * 	\brief Print the warnings from the buffer or do thing if the buffer is empty.
	 */
	~Warning();

	/**
	 * \brief Initialize the warning and the internal buffer.
	 *
	 * @param message The message the internal buffer is set to.
	 */
	Warning(std::string const& message);

	/**
	 * \brief Allow the user to insert any type into the warning.
	 *
	 * It can be used with any type that can be stringified using the operator<< of a stringstream
	 *
	 * @param message The message the internal buffer is set to.
	 */
	template<typename T>
	Warning& operator<< (T const& message);

	void print() const;
private:

	std::string _buffer;
	std::stringstream _sstrBuff;
};

} /* namespace error_handling */
} /* namespace scallop */
#include "scallop/error_handling/src/Warning.hpp"
#endif /* SCALLOP_ERROR_HANDLING_WARNING_H_ */
