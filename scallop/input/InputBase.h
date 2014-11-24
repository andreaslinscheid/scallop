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

	/**
	 * Write the internally stored input manual to the fileName.
	 *
	 * @param fileName The full path of the manual to be created.
	 */
	void build_input_manual(std::string const& fileName) const;

	/**
	 * Take an input file and fill the internal variables with content.
	 *
	 * @param inputFile The parsed input file.
	 */
	void parse_variables(InputFile const& inputFile);

protected:

	bool _isInit = false;

	template<typename T>
	void get_option(std::string const& valueString,T const& defaultValue,T &value) const;

	template<typename T>
	void get_option(std::string const& valueString,T &value) const;

	template<typename T>
	void get_option(std::string const& valueString,std::vector<T> const& defaultValue,std::vector<T> &values) const;

	template<typename T>
	void get_option(std::string const& valueString,std::vector<T> &values) const;

	void add_option_to_manual(std::string const& textWithOptionDescription);

	void register_variable_parsing( void( derived:: * function )(InputFile const&) );

private:

	std::string _manual;

	std::vector< void(derived::*)(InputFile const&) > _mebrFctnPtrToVariableParsing;
};

/**
 * Allows to add comma in macro expansion values.
 *
 * Use for example {0 COMMA_SUBSTITUTION 1 COMMA_SUBSTITUTION 2 COMMA_SUBSTITUTION 3}
 * to initialize, say, a vector with the list {0,1,2,3} via a macro variable
 */
#define COMMA_SUBSTITUTION ,

/**
 * This macro allows to easily add new input variables with a default value.
 *
 *
 * To add a new input variable use
 * 		INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
 * 			quantityName,
 * 			description,
 * 			defaultValueText,
 * 			defaultValueStatement,
 * 			typeQ);
 *
 *	Here:
 *\param quantityName  		Is the name of quantity in the input file and
 *							internally in the Input module
 *\param description		A description of the variable that appears in the manual.
 *\param defaultValueText	A description text of the default value such as "zero" or "0" or "sequence 0 to 3".
 *\param defaultValueStatement	A compiling expression for the default such as 0 or
 *							{0 COMMA_SUBSTITUTION 1 COMMA_SUBSTITUTION 2 COMMA_SUBSTITUTION 3}
 *\param typeQ				The type of the variable such as double or size_t or std::vector<size_t>
 *
 */
#define INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(											\
		quantityName,description,defaultValueText,defaultValueStatement,typeQ)				\
public:                                                                                 	\
	typeQ const& get_##quantityName() const {                                            	\
		if ( not this->_isInit ){															\
			error_handling::Error("Calling get_"#quantityName"()"							\
									"before calling parse_variables()",1);					\
		};																					\
		return _##quantityName;                                                         	\
	};                                                                                  	\
private:																					\
	typeQ _##quantityName;																	\
	void parse_variable_##quantityName(InputFile const& inputFile) {						\
		typeQ defaultValue = defaultValueStatement;											\
		std::string valueInInputFile = inputFile.get_input_config_value(#quantityName);		\
		this->get_option(valueInInputFile,defaultValue, _##quantityName);					\
	}																						\
	void * perform_##quantityName##_processing() {											\
		this->add_option_to_manual( std::string()+											\
		  "\n========================================================================\n"	\
			"||Variable :             || "+#quantityName+"\n"								\
			"||-----------------------||---------------------------------------------\n"	\
			"||Type is :              || "+#typeQ+"\n"										\
			"||-----------------------||---------------------------------------------\n"	\
			"||The default value is : || "+defaultValueText+"\n"							\
			"========================================================================\n"	\
			"|Description:|\n"																\
			"--------------\n"																\
			""+description+"\n"                                      						\
			"========================================================================\n\n");\
		this->register_variable_parsing(													\
			&std::remove_pointer<decltype(this)>::type ::parse_variable_##quantityName);	\
		return 0;																			\
	};																						\
	void * call_##quantityName##_processing													\
				= this->perform_##quantityName##_processing()								\

/**
 * 	Similar macro as \ref INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT without default values.
 *
 * 	This requires the variable to be set in the input process or the program errors out.
 * 	See documentation of \ref INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT
 */
#define INPUTBASE_INPUT_OPTION_MACRO(quantityName,description,typeQ)						\
public:                                                                                 	\
	typeQ const& get_##quantityName() const {                                            	\
		if ( not this->_isInit ){															\
			error_handling::Error("Calling get_"#quantityName"()"							\
									"before calling parse_variables()",1);					\
		};																					\
		return _##quantityName;                                                         	\
	};                                                                                  	\
private:																					\
	typeQ _##quantityName;																	\
	void parse_variable_##quantityName(InputFile const& inputFile) {						\
		std::string valueInInputFile = inputFile.get_input_config_value(#quantityName);		\
		this->get_option(valueInInputFile,_##quantityName);									\
	}																						\
	void * perform_##quantityName##_processing() {											\
		this->add_option_to_manual( std::string()+											\
		  "\n========================================================================\n"	\
			"||Variable :             || "+#quantityName+"\n"								\
			"||-----------------------||---------------------------------------------\n"	\
			"||Type is :              || "+#typeQ+"\n"										\
			"||-----------------------||---------------------------------------------\n"	\
			"||The default value :    || This variable is MANDATORY!\n"						\
			"========================================================================\n"	\
			"|Description:|\n"																\
			"--------------\n"																\
			""+description+"\n"                                      						\
			"========================================================================\n\n");\
		this->register_variable_parsing(													\
			&std::remove_pointer<decltype(this)>::type ::parse_variable_##quantityName);	\
		return 0;																			\
	};																						\
	void * call_##quantityName##_processing													\
				= this->perform_##quantityName##_processing()								\


} /* namespace input */
} /* namespace scallop */
#include "scallop/input/src/InputBase.hpp"
#endif /* SCALLOP_INPUT_INPUTBASE_H_ */
