/*	This file Error.cpp is part of scallop.
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

#include "scallop/error_handling/Error.h"
#include <string>
#include <cstdlib>
#include <iostream>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>

namespace scallop {
namespace error_handling {

Error::Error(int errorCode){
	this->call_error_no_ref("",errorCode);
}

Error::Error(std::string const& description,int errorCode){
	this->call_error_no_ref(description,errorCode);
}

Error::Error(int line,std::string const& file,std::string const& description,size_t errorCode) {
	this->call_error_with_code_ref(line,file,description,static_cast<int>(errorCode));
}

void Error::call_error_with_code_ref(int line,
		std::string const& file,
		std::string const& description,
		int errorCode){

	std::cout << "In line " << line << " of file " << file << " :" << std::endl;
	this->call_error_no_ref(description,errorCode);
}
void Error::call_error_no_ref(std::string const& description,int errorCode){
	if (errorCode == 0)
		return;
	std::cout << "ERROR occurred!";
	if ( not description.empty() )
		std::cout << "\tProblem description : " <<description << std::endl;

	//If we are using debugging symbols generate a stack trace
#ifdef DEBUG_BUILD
	handler(errorCode);
#else
	std::exit(errorCode);
#endif
}

//backtrace and backtrace_symbols_fd do not malloc and should thus be ok to call on throw.
void signal_handler(int signal) {
	 void * stackEntries[20];
	  size_t size;

	  //get pointers to all entries on the stack
	  size = backtrace(stackEntries, 20);

	  // print out all the frames to stderr
	  std::cerr << "Error: signal %d:\n" << signal << std::endl;
	  backtrace_symbols_fd(stackEntries, size, STDERR_FILENO);
	  std::exit(1);
}

void Error::handler(int signal) {
	signal_handler(signal);
}

//install signal handler
void Error::generate_stacktrace_on_segfault() {
#ifdef DEBUG_BUILD
	signal(SIGSEGV, signal_handler);
#endif
}

} /* namespace error_handling */
} /* namespace scallop */
