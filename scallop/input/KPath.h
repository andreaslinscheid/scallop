/*	This file KPath.h is part of scallop.
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
 *  Created on: Nov 27, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_INPUT_KPATH_H_
#define SCALLOP_INPUT_KPATH_H_

namespace scallop
{
namespace input
{

template<typename T>
class KPath
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;

	typedef typename auxillary::TemplateTypedefs<bT>::scallop_vector V;

	void read_kpath_file( std::string const& filename );

	V const& get_k_path() const;

	template<class DT, class OT, typename T2>
	void band_structure_gnuplot(
			DT const& data,
			OT const& omega,
			parallel::IrregularGridDistribution<T2> const& irrGrid,
			std::string const& filename,
			std::string quantityLabel) const;
private:

	size_t dim_ = 0;

	V kptsTotal_;

	//Note that we don't need labels_ on processors other than the ioproc
	std::vector<std::pair<size_t,std::string> > labels_;
};

} /* namespace input */
} /* namespace scallop */

#include "scallop/input/src/KPath.hpp"
#endif /* SCALLOP_INPUT_KPATH_H_ */
