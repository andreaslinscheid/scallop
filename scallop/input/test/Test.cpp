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

#include "scallop/input/test/Test.h"
#include "scallop/input/input.h"
#include "scallop/output/TextFile.h"
#include "scallop/input/test/InputBaseTest.h"
#include "scallop/output/TerminalOut.h"

namespace scallop {
namespace input {
namespace test {

void Test::run_test() {

	test_input_file_parsing();
}

void Test::test_input_file_parsing() {

	output::TerminalOut msg;
	std::string tmpFolder = "/tmp/";

	//create an example input file
	std::string const testInputFileName = tmpFolder+"testin.dat";
	create_test_input_file(testInputFileName);

	//do as if this file was given as a option during program start.
	char ** argv;
	argv = new char* [2];
	argv[1] = new char [testInputFileName.size()+1];
	std::copy(testInputFileName.c_str(),testInputFileName.c_str()+testInputFileName.size(),argv[1]);
	argv[1][testInputFileName.size()] = '\0';

	//test the startup routines and parse the test input file
	Setup setup(2,argv);

	//get the options from the input file
	InputBaseTest inputFileTest;
	inputFileTest.parse_variables( setup.get_parsed_input_file() );

	//Test the manual generation
	std::string const testManualName = tmpFolder + "manualTest.txt";
	inputFileTest.build_input_manual(testManualName);

	//Print unused parameters
	std::vector<std::string> unreadOptions =
			setup.get_parsed_input_file().get_list_unread_input_parameters();
	msg << "Testing that the following input parameters were not used: ";
	for ( auto &opt  : unreadOptions )
		msg << '\t'<< opt;

	msg << "Read in sizeTest : " << inputFileTest.get_sizeTest();
	msg << "Read in doubleTest : " << inputFileTest.get_doubleTest();
	msg << "Read in boolTest : " << std::boolalpha << inputFileTest.get_boolTest();
	msg << "Read in vectorSizeT with the elements : ";
	std::vector<size_t> vect = inputFileTest.get_vectorSizeT();
	for ( auto & e : vect )
		msg << e;
}

void Test::create_test_input_file(std::string const& fileName){
	std::string const inputFileContent = "#Testfile with a comment as header line\n"
			"#\t a tab and some test input data:\n"
			"sizeTest=42\ndoubleTest=3.4\n\n"
			"vectorSizeT=1 2 3 !testing vector reading and an inline comment\n"
			"#to redefine the same key but only with the same value.\n"
			"#In the next line we test an empty line:\n\n"
			"boolThatWillNotBeUsed = true #The program should inform about this unused variable!\n"
			"boolTest = false";

	output::TextFile textFile;
	textFile.write(fileName,inputFileContent);
}

} /* namespace test */
} /* namespace input */
} /* namespace scallop */
