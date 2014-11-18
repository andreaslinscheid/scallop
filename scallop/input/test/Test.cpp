/*	This file Test.cpp is part of scallop.
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
 *  Created on: Nov 15, 2014
 *      Author: Andreas Linscheid
 */

#include "Test.h"
#include "scallop/input/input.h"
#include "scallop/output/TextFile.h"
#include "scallop/input/test/InputBaseTest.h"

namespace scallop {
namespace input {
namespace test {

void Test::run_test() {

	test_input_file_paring();
}

void Test::test_input_file_paring() {

	//create an example input file
	std::string const testInputFileName = "testin.dat";
	create_test_input_file(testInputFileName);

	//do as if this file was given as a option during program start.
	char ** argv;
	argv = new char* [2];
	argv[1] = new char [testInputFileName.size()];
	std::copy(testInputFileName.c_str(),testInputFileName.c_str()+testInputFileName.size(),argv[1]);

	//test the startup routines and parse the test input file
	Setup setup(2,argv);

	//get the options from the input file
	InputBaseTest inputFileTest;
	inputFileTest.parse_all_registered_variables( setup.get_parsed_input_file() );

	std::string const testManualName = "manualTest.txt";
	inputFileTest.build_input_manual(testManualName);
}

void Test::create_test_input_file(std::string const& fileName){
	std::string const inputFileContent = "#Testfile with a comment as header line\n"
			"#\t a tab and some test input data:"
			"sizeTest=42\ndoubleTest=3.4\n";

	output::TextFile textFile;
	textFile.write(fileName,inputFileContent);
}

} /* namespace test */
} /* namespace input */
} /* namespace scallop */
