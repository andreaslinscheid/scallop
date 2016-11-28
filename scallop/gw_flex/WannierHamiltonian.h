/*	This file WannierHamiltonian.h is part of scallop.
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
 *  Created on: Nov 26, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_WANNIERHAMILTONIAN_H_
#define SCALLOP_GW_FLEX_WANNIERHAMILTONIAN_H_

namespace scallop
{
namespace gw_flex
{

template<typename T>
class WannierHamiltonian
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;
	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector v;
	typedef typename auxillary::TemplateTypedefs<bT>::scallop_vector vbT;

	void load_wan_ham( std::string const & fileNameWannierHam );

	void compute_at_k( vbT kpts, size_t nkpts, v & unitary, vbT & energyEV ) const;

	size_t get_nOrb() const;

private:

	size_t nOrb_ = 0;

	vbT RGrid_;

	vbT weights_;

	v wanHam_;

};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/WannierHamiltonian.hpp"
#endif /* SCALLOP_GW_FLEX_WANNIERHAMILTONIAN_H_ */
