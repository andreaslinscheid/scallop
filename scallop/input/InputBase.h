/*	This file InputBase.h is part of scallop.
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

#ifndef SCALLOP_INPUT_INPUTBASE_H_
#define SCALLOP_INPUT_INPUTBASE_H_

#include "scallop/input/InputFile.h"
#include <string>
#include <vector>
#include <type_traits>

namespace scallop {
namespace input {

template<class derived>
class InputBase {

public:

	template<typename T>
	void get_option(std::string const& valueString,T const& defaultValue,T &value) const;

	template<typename T>
	void get_option(std::string const& valueString,T &value) const;

	template<typename T>
	void get_option(std::string const& valueString,std::vector<T> const& defaultValue,std::vector<T> &values) const;

	template<typename T>
	void get_option(std::string const& valueString,std::vector<T> &values) const;

//	void get_option(std::string const& valueString,
//			std::vector<bool> const& defaultValue,
//			std::vector<bool> &values) const ;

	void build_input_manual(std::string const& fileName) const;

	void add_option_to_manual(std::string const& textWithOptionDescription);

	void parse_all_registered_variables(InputFile const& inputFile);

	void register_variable_parsing( void(derived::* function )(InputFile const&) );
private:

	std::vector< void(derived::*)(InputFile const&) > _mebrFctnPtrToVariableParsing;

	std::string _manual;
};

/**
 * This macro allows to easily add new input variables.
 *
 * To add a new input variable
 *
 */
#define INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(\
		quantityName,description,defaultValueText,defaultValueStatement,typeQ)\
public:                                                                                 	\
	typeQ const& get_##quantityName() const {                                            	\
		return _##quantityName;                                                         	\
	};                                                                                  	\
private:																					\
	typeQ _##quantityName;																	\
	void parse_variable_##quantityName(InputFile const& inputFile) {													\
		typeQ defaultValue = defaultValueStatement;											\
		std::string valueInInputFile = inputFile.get_input_config_value(#quantityName);		\
		this->get_option(valueInInputFile,defaultValue, _##quantityName);					\
	}																						\
	void * perform_##quantityName##_processing() {											\
		this->add_option_to_manual( std::string()+											\
				"\n==========================================\n"							\
				"Variable "+#quantityName+" of type "+#typeQ									\
					+"\nThe default is "+#defaultValueText+"\n"								\
				"Description:\n"+#description+" \n"                                      	\
				"==========================================\n\n");                     		\
		this->register_variable_parsing(&std::remove_pointer<decltype(this)>::type ::parse_variable_##quantityName);					\
		return 0;																			\
	};																						\
	void * call_##quantityName##_processing													\
				= this->perform_##quantityName##_processing()								\

#define INPUTBASE_INPUT_OPTION_MACRO(quantityName,description,type)							\
public:                                                                                 	\
	type const& get_##quantityName() const {                                            	\
		return _##quantityName;                                                         	\
	};                                                                                  	\
private:																					\
	type _##quantityName;																	\
	void * add_##quantityName##_description() {												\
			this->add_option_to_manual( std::string()+										\
				"\n==========================================\n"							\
				"Variable "+#quantityName+" of type "+#type									\
					+"\nThis variable is REQUIRED on input!\n"								\
				"Description:\n"+#description+" \n"                                      	\
				"==========================================\n\n");                     		\
		return 0;																			\
	};																						\
	void * call_##quantityName##_description_add											\
											= add_##quantityName##_description()			\

} /* namespace input */
} /* namespace scallop */
#include "scallop/input/src/InputBase.hpp"
#endif /* SCALLOP_INPUT_INPUTBASE_H_ */
