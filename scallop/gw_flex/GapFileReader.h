/*	This file GapFileReader.h is part of scallop.
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
 *  Created on: Jan 18, 2017
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_GAPFILEREADER_H_
#define SCALLOP_GW_FLEX_GAPFILEREADER_H_

#include "scallop/auxillary/TypeMapComplex.h"
#include "scallop/error_handling/Error.h"
#include <string>
#include <vector>

namespace scallop
{
namespace gw_flex
{

class GapFileReader
{
public:

	void read_file(std::string const& filename );

	template<typename T>
	std::complex<T> operator() (std::vector<T> const& k,
			T omegan,
			size_t component,
			size_t orbitalIndex) const;
private:

	typedef double (GapFileReader::* Shape_t )( std::vector<double> const&) const;

	typedef struct
	{
		size_t spinState_ = 0;
		double width_ = 0;
		std::vector<double> orbFactors_;
		Shape_t shape_ = NULL;

	} ChannelFuncPara;

	std::vector<ChannelFuncPara> channels_;

	double s_wave( std::vector<double> const& k ) const;

	double px_wave( std::vector<double> const& k ) const;
	double py_wave( std::vector<double> const& k ) const;
	double pz_wave( std::vector<double> const& k ) const;

	double dxy_wave( std::vector<double> const& k ) const;
	double dyz_wave( std::vector<double> const& k ) const;
	double dz2_wave( std::vector<double> const& k ) const;
	double dxz_wave( std::vector<double> const& k ) const;
	double dx2my2_wave( std::vector<double> const& k ) const;

	double radius( std::vector<double> const& k ) const;
	double radius2( std::vector<double> const& k ) const;
};


template<typename T>
std::complex<T> GapFileReader::operator() (
		std::vector< T > const& k,
		T omegan,
		size_t component,
		size_t orbitalIndex) const
{
	#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

	std::vector<double> kFirstZone( k.begin(), k.end() );
	for (auto &kfz : kFirstZone )
	{
		kfz -= std::floor(kfz+0.5);
		kfz *= M_PI;
	}

	std::complex<T> result = std::complex<T>(0);
	for ( auto c : channels_ )
	{
		if ( ! (component == c.spinState_ ) )
			continue;
		if ( orbitalIndex >= c.orbFactors_.size() )
			error_handling::Error("Requesting gap initialization for values not specified on input");
		double g = CALL_MEMBER_FN(*this,(c.shape_))(kFirstZone)*c.orbFactors_[orbitalIndex];
		g *= c.width_*c.width_/(omegan*omegan+c.width_*c.width_);
		result += std::complex<T>(g);
	}
	return result;
}
} /* namespace gw_flex */
} /* namespace scallop */

#endif /* SCALLOP_GW_FLEX_GAPFILEREADER_H_ */
