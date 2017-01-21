/*	This file Driver.h is part of scallop.
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
 *  Created on: Dec 13, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_DRIVER_H_
#define SCALLOP_GW_FLEX_DRIVER_H_

#include "scallop/input/Configuration.h"
#include "scallop/gw_flex/KohnShamBandStructure.h"
#include "scallop/gw_flex/KohnShamGreensFunctionOrbital.h"
#include "scallop/gw_flex/SelfEnergy.h"
#include "scallop/gw_flex/ChemicalPotentialShifting.h"
#include "scallop/gw_flex/InteractionMatrix.h"
#include "scallop/gw_flex/SpinSusceptibility.h"
#include "scallop/gw_flex/ChargeSusceptibility.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
class Driver
{
public:

	Driver(input::Configuration config);

	void converge();

private:

	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector V;

	typedef typename auxillary::TypeMapComplex<T>::type bT;

	input::Configuration config_;

	size_t iter_ = 0;

	KohnShamBandStructure<T> elstr_;

	V matsFreq_;
	V prevMatsFreq_;

	InteractionMatrix<T> Isf_;

	InteractionMatrix<T> Ic_;

	KohnShamGreensFunctionOrbital<T> ksGF_;

	KohnShamGreensFunctionOrbital<T> ksGT_;

	GreensFunctionOrbital<T> G_;

	ChemicalPotentialShifting<T> chmpot_;

	SpinSusceptibility<T> suscSF_;

	ChargeSusceptibility<T> suscC_;

	SelfEnergy<T> se_;

	SelfEnergy<T> seL_;

	typename GeneralizedSusceptibility<T>::AdiabaticUpscale spinAdiabaticScale_;

	void initialize_this_T_N( double temp, double & Ne );

	void set_gap_symmetry_breaking( std::string const& filename );

	void prepare_iteration( double temp, bool set_GF_from_KS, bool set_GF_from_SE );

	void shift_chemical_pot( double temp, double Ne );

	void compute_eff_interactions( bT beta );

	void compute_self_energies();

	bool check_convergence() const;

	void mix_iterations();

	void solve_Dyson();

	void report( bT beta ) const;

	void post_process();
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/Driver.hpp"
#endif /* SCALLOP_GW_FLEX_DRIVER_H_ */
