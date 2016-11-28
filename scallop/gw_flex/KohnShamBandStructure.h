/*	This file KohnShamBandStructure.h is part of scallop.
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
 *  Created on: Nov 25, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_KOHNSHAMBANDSTRUCTURE_H_
#define SCALLOP_GW_FLEX_KOHNSHAMBANDSTRUCTURE_H_

#include "scallop/gw_flex/UnitaryWannierKSBands.h"
#include "scallop/gw_flex/WannierHamiltonian.h"
#include "scallop/gw_flex/BandstructureModels.h"

namespace scallop
{
namespace gw_flex
{

template<typename T >
class KohnShamBandStructure
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;
	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector v;
	typedef typename auxillary::TemplateTypedefs<bT>::scallop_vector vbt;

	void initialize_from_file(
			std::vector<size_t> grid,
			std::string const & fileName );

	bT operator() (size_t ik, size_t n) const;

	UnitaryWannierKSBands<T> const & get_unitary() const;

	vbt const & get_bands() const;

	void compute_at_k( vbt const& kpoints,size_t nK,vbt & enk_,v unitary) const;

	size_t get_nOrb() const;

	parallel::GridDistribution<T> const& get_spaceGrid_proc() const;
private:

	bool useModel_ = false;

	vbt enk_;

	UnitaryWannierKSBands<T> akil_;

	WannierHamiltonian<T> wanHam_;

	std::shared_ptr<BandstructureModels<T> > model_ = NULL;

	parallel::GridDistribution<T> gridDistr_;

	void set_model(std::istream & stream);
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/KohnShamBandStructure.hpp"
#endif /* SCALLOP_GW_FLEX_KOHNSHAMBANDSTRUCTURE_H_ */
