/*	This file BasicFunctions.cpp is part of scallop.
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
 *  Created on: Dec 8, 2016
 *      Author: A. Linscheid
 */

#include "scallop/auxillary/BasicFunctions.h"
#include <string>
#include <sys/types.h>
#include <unistd.h>
#include <regex>

namespace scallop
{
namespace auxillary
{

std::string BasicFunctions::get_scallop_path() const
{
	std::string fullFileName = "";

	// Code taken from: http://www.gamedev.net/community/forums/topic.asp?topic_id=459511
	std::string pathToScallop = "";
	pid_t pid = getpid();
	char buf[20] = {0};
	sprintf(buf,"%d",pid);
	std::string _link = "/proc/";
	_link.append( buf );
	_link.append( "/exe");
	char proc[512];
	int ch = readlink(_link.c_str(),proc,512);
	if (ch != -1) {
		proc[ch] = 0;
		pathToScallop = proc;
		std::string::size_type t = pathToScallop.find_last_of("/");
		pathToScallop = pathToScallop.substr(0,t);
	}

	pathToScallop = std::regex_replace(pathToScallop, std::regex("(/tests)?/bin\\..*/?"), "");

	return fullFileName = pathToScallop + std::string("/");
}

std::complex<double> BasicFunctions::singlet_Re_channel(size_t a1, size_t s1, size_t a2, size_t s2)
{
	return -pauli_y(a1,a2)*pauli_y(s1,s2);
}

std::complex<double> BasicFunctions::singlet_Im_channel(size_t a1, size_t s1, size_t a2, size_t s2)
{
	return pauli_x(a1,a2)*std::complex<double>(0,1)*pauli_y(s1,s2);
}

std::complex<double> BasicFunctions::triplet_x_Re_channel(size_t a1, size_t s1, size_t a2, size_t s2)
{
	return -std::complex<double>(0,1)*pauli_y(a1,a2)*pauli_z(s1,s2);
}

std::complex<double> BasicFunctions::triplet_x_Im_channel(size_t a1, size_t s1, size_t a2, size_t s2)
{
	return -pauli_x(a1,a2)*pauli_z(s1,s2);
}

std::complex<double> BasicFunctions::triplet_y_Re_channel(size_t a1, size_t s1, size_t a2, size_t s2)
{
	return -std::complex<double>(0,1)*pauli_y(a1,a2)*pauli_0(s1,s2);
}

std::complex<double> BasicFunctions::triplet_y_Im_channel(size_t a1, size_t s1, size_t a2, size_t s2)
{
	return -pauli_x(a1,a2)*pauli_0(s1,s2);
}

std::complex<double> BasicFunctions::triplet_z_Re_channel(size_t a1, size_t s1, size_t a2, size_t s2)
{
	return std::complex<double>(0,1)*pauli_y(a1,a2)*pauli_x(s1,s2);
}

std::complex<double> BasicFunctions::triplet_z_Im_channel(size_t a1, size_t s1, size_t a2, size_t s2)
{
	return pauli_x(a1,a2)*pauli_x(s1,s2);
}

std::complex<double> BasicFunctions::pauli_0(size_t i, size_t j)
{
	assert( (i < 2) && (j < 2) );
	return i == j ? std::complex<double>(1.0) : std::complex<double>(0.0);
}


std::complex<double> BasicFunctions::pauli_x(size_t i, size_t j)
{
	assert( (i < 2) && (j < 2) );
	return i != j ? std::complex<double>(1.0) : std::complex<double>(0.0);
}

std::complex<double> BasicFunctions::pauli_y(size_t i, size_t j)
{
	assert( (i < 2) && (j < 2) );
	return i != j ? std::complex<double>(0.0,1.0)*(j==1?-1.0:1.0) : std::complex<double>(0.0);
}

std::complex<double> BasicFunctions::pauli_z(size_t i, size_t j)
{
	assert( (i < 2) && (j < 2) );
	return i == j ? std::complex<double>(1.0)*(i==1?-1.0:1.0) : std::complex<double>(0.0);
}

} /* namespace auxillary */
} /* namespace scallop */
