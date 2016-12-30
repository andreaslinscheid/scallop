/*	This file ChemicalPotentialShifting.h is part of scallop.
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
 *  Created on: Dec 15, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_CHEMICALPOTENTIALSHIFTING_H_
#define SCALLOP_GW_FLEX_CHEMICALPOTENTIALSHIFTING_H_

#include "scallop/auxillary/TemplateTypedefs.h"
#include "scallop/gw_flex/KohnShamBandStructure.h"
#include "scallop/gw_flex/KohnShamGreensFunctionOrbital.h"
#include "scallop/gw_flex/GreensFunctionOrbital.h"
#include "scallop/auxillary/LinearAlgebraInterface.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
class ChemicalPotentialShifting
{
public:

	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector V;

	typedef typename auxillary::TypeMapComplex<T>::type bT;

	void determine_new_chemPot(
			bT Nelectrons,
			bT invTemp,
			KohnShamBandStructure<T> const& elstr,
			KohnShamGreensFunctionOrbital<T> const& gks,
			GreensFunctionOrbital<T> const& g,
			std::pair<bool,bT> & newChemPot);

private:

	V buffer1_,buffer2_;

	void compute_N_elec_mB(
			KohnShamGreensFunctionOrbital<T> const& gks,
			GreensFunctionOrbital<T> const& g,
			bT beta,
			bT muDelta,
			bT & nElec,
			MemoryLayout & gfl);
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/ChemicalPotentialShifting.hpp"
#endif /* SCALLOP_GW_FLEX_CHEMICALPOTENTIALSHIFTING_H_ */
