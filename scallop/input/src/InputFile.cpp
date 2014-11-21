/*	This file InputFile.cpp is part of scallop.
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

#include "scallop/input/InputFile.h"
#include "scallop/error_handling/Error.h"
#include <fstream>
#include <vector>
#include <string>

namespace scallop {
namespace input {

void InputFile::read_input_file(
		std::string const & fileName,
		std::string & infileContent) const {
	infileContent.clear();

	//open the file
	std::ifstream file;
	file.open(fileName.c_str(),std::ios::in);
	if ( not file.good() ){
		std::string const msg = std::string("File ") + fileName + std::string(" not readable!");
		error_handling::Error(msg, 1);
	}

	//read the file
	std::string line;
	while ( std::getline(file,line) ) {
		//Check if the last char is \r as might occur with vim edited files
		if (line.c_str()[line.size()-1] == '\r') {
			line.resize(line.size () - 1);
		}
		infileContent.append(line);
		infileContent.append("\n");
		line.clear();
	}
	file.close();
}


//The strategy is to break the input into parts separated by '='.
// except the first and the last element this means we have a block with
//	value+key. Keywords are not separated by whitespace.
void InputFile::parse_input(std::string const & input) {
	//remove comments

	//break input into key/value pairs
    std::stringstream ss(input);
    std::string item;
    std::vector<std::string> elements;
    while( std::getline(ss, item, '=') ) {
        elements.push_back(item);
    }

    std::vector<std::pair<std::string,std::string> > keyValuePairs;
    auto it = elements.begin();
    std::string value,key = *it;
    for ( ++it; it != (--elements.end()); ++it ) {

    	std::vector<std::string> data;
        std::stringstream valuePlusNextKeyRaw(*it);
        while( std::getline(valuePlusNextKeyRaw, item) ) {
        	data.push_back(item);
        }

        value.clear();
        for ( std::vector<std::string>::iterator itd = data.begin(); itd != (-- data.end()); ++itd )
        	value += *itd + ' ';

        trim_string(key);
        trim_string(value);
    	keyValuePairs.push_back( std::make_pair(key,value) );

    	key = data.back();
    }
	keyValuePairs.push_back( std::make_pair(key,elements.back()) );

	//insert the key and value pairs into the map.
	for ( auto const& elmentPair : keyValuePairs ) {

		auto ret = _inputFileKeyValue.insert(elmentPair);
		if ( not ret.second ) { //i.e. the key exists already

			if (  (*ret.first).second.compare(elmentPair.second) != 0 ) {
				std::string const msg = std::string("Key ") + (*ret.first).first
						+ " appears twice with disagreeing value, i.e. first " + (*ret.first).second
						+ " and second with " + elmentPair.second;
				error_handling::Error(msg , 1);
			}

		}

		std::pair<std::string,bool> keyReadThisValue;
		keyReadThisValue.first = elmentPair.first;
		keyReadThisValue.second = false; //has not yet been read
		_keyWasRead.insert(keyReadThisValue);
	}
}

std::string InputFile::get_input_config_value( std::string const& key ) const {
	std::string result;
	auto it = _inputFileKeyValue.find(key);
	if ( it != _inputFileKeyValue.end() ){
		result = (*it).second;

		//keeping track of the read input options. Note that _keyWasRead is mutable
		_keyWasRead[key] = true;
	}

	return result;
}


void InputFile::trim_string(std::string & str,
		std::string const& whitespace) const {
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
    	str = ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    std::string tempStr = str.substr(strBegin, strRange);
    std::swap(tempStr,str);
}


std::vector<std::string> InputFile::get_list_unread_input_parameters() const {
	std::vector<std::string> result;
	for ( auto &element : _keyWasRead ) {
		if ( not element.second )
			result.push_back(element.first);
	}
	return result;
}

} /* namespace input */
} /* namespace scallop */
