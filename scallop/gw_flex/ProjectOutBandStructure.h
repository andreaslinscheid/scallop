/*	This file ProjectOutBandStructure.h is part of scallop.
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

#ifndef SCALLOP_GW_FLEX_PROJECTOUTBANDSTRUCTURE_H_
#define SCALLOP_GW_FLEX_PROJECTOUTBANDSTRUCTURE_H_

#include "scallop/gw_flex/KohnShamBandStructure.h"
#include "scallop/gw_flex/SelfEnergy.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
class ProjectOutBandStructure
{
public:
	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector V;
	typedef typename auxillary::TemplateTypedefs< std::complex<float> >::scallop_vector container_type;

	ProjectOutBandStructure(
			KohnShamBandStructure<T> bands);

	size_t size_per_block() const;

	size_t size_per_block_after() const;

	void project_out_KS_bands_half_freq(
			SelfEnergy<T> const& SE,
			V const& MatsFreq,
			V & mappedMatsFreq, container_type & mappedData) const;

	void before( container_type & selfEnergyAlongPath,
			parallel::IrregularGridDistribution<T> const& path,
			V const& MatsGrid) const;

	void after( container_type & dataPath,
			parallel::IrregularGridDistribution<T> const& path,
			V const& omega) const;
private:

	KohnShamBandStructure<T> bands_;
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/ProjectOutBandStructure.hpp"
#endif /* SCALLOP_GW_FLEX_PROJECTOUTBANDSTRUCTURE_H_ */
