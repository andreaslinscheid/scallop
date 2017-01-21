/*	This file Configuration.h is part of scallop.
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
 *  Created on: Nov 25, 2014
 *      Author: Andreas Linscheid
 */

#ifndef SCALLOP_INPUT_CONFIGURATION_H_
#define SCALLOP_INPUT_CONFIGURATION_H_

#include "scallop/input/InputBase.h"

namespace scallop {
namespace input {

class Configuration : public InputBase<Configuration>
{
public:

	INPUTBASE_INPUT_OPTION_MACRO(
			f_elstr,
			"The file that determines the electronic structure.\n"
			"Could be a Hamiltonian file in the wannier90 format or a model file.\n"
			"Must be set.\n",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			grid,
			"The grid used for the FLEX algorithm, i.e. self-energy etc.\n"
			"Can be 1-3D with parallelization supported for 2 and 3D.\n"
			"The number of points must be > number of processors.\n",
			"10x10",
			{10 COMMA_SUBSTITUTION 10},
			std::vector<size_t>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			method,
			"Decide which method is used to compute SC\n"
			"Possible choices are:\n\t isotropic Eliashberg 'eli'\n"
			"\t orbital resolved Gw-FLEX 'kgw'",
			"kgw",
			"kgw",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			adj_mu,
			"Re-adjust the chemical potential such that the number of electrons in the system is fix\n",
			"true",
			true,
			bool);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			maxIter,
			"Maximal number of iterations\n",
			"1000",
			1000,
			size_t);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			temp,
			"The temperatures where to run\n"
			"in Kelvin.",
			"10.0",
			{10.0},
			std::vector<double>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			nelec,
			"The number of electrons per spin in the system\n"
			"The chemical potential will be determined on startup, such that this numbers are reached.\n"
			"A negative number indicates that the chemical potential from the input electronic structure is used.",
			"-1.0",
			{-1.0},
			std::vector<double>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			MCut,
			"Cutoff for the Matsubara frequencies in the self-energy\n"
			"A negative number indicates that we choose the band width",
			"-1.0",
			-1.0,
			double);


	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			f_sf,
			"The file with the interaction matrix elements for the spin fluctuations.\n"
			"The default is empty, meaning no spin fluctuations will be considered.\n",
			"",
			"",
			std::string);


	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			f_c,
			"The file with the interaction matrix elements for the charge fluctuations.\n"
			"The default is empty, meaning no charge fluctuations will be considered.\n",
			"",
			"",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			acc_d,
			"Convergence threshold for self-energies. The calculation is considered converged if all data points of the self-energy\n"
			"differ by no more than these number of significant digits.",
			"3",
			3,
			int);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			mix,
			"mixing parameter for linear mixing of Self-energies.\n"
			"Range 0-1 where 1 is 100% of the old iteration.",
			"0.6",
			0.6,
			double);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			cmp_zero,
			"Value below which a self-energy part is considered zero\n"
			"Unit is meV",
			"1e-5",
			1e-5,
			double);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			vbl,
			"Verbosity level to indicate what detail should be written to STDOUT\n"
			"Possible values are 'l' (low), 'm' (medium) or 'h' (high)",
			"m",
			"m",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			f_spec,
			"The file name template of the single particle spectral function.\n"
			"Setting this file means that the spectral function will be computed and written to the file '<f_spec>.dat'.\n"
			"The code will also generate a gnuplot script under '<f_spec>.gp' to plot the data.\n"
			"If this option is set, f_kpath must be set too.\n",
			"",
			"",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			r_omega,
			"The nearly real axis frequency grid.\n"
			"The code expects four values and will abort if a different number is set.\n"
			"The first (second) is the lower (upper) bound in meV and the third is \n"
			"converted to integer for the number of points in the grid while the fourth is the imaginary offset from the real axis.\n",
			"The band edges +10% bandwidth of the non-interacting system, 1000 pts, and 1e-5*bandwidth.",
			{},
			std::vector<double>);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			f_kpath,
			"The file name with a k path that is used by postprocessing things.\n",
			"",
			"",
			std::string);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			maxDownscale,
			"If the code hits an instability in the RPA enhancement of the susceptiblity, we apply a downscaling \n"
			"of all channels of the respective interaction up to this value. The code tries (scale+maxDownscale)/2"
			"starting from scale=1 'scale_res' times until we are in the stabel regime.",
			"0.25",
			0.25,
			double);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			scale_res,
			"see 'maxDownscale'",
			"10",
			10,
			size_t);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			scale_SE_iter,
			"Once a scaling succeeds, i.e. we have obtained a stable interaction, we compute the self-energy."
			"The code will try 'scale_SE_iter' times to reduce the scaling back to 1 before giving up.",
			"10",
			10,
			size_t);

	INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
			f_gap_i,
			"The file with the starting gap.\n"
			"If the file is not set, we assume a vanishing SC phase, i.e. the symmetry will not be broken",
			"",
			"",
			std::string);

};

} /* namespace input */
} /* namespace scallop */
#endif /* SCALLOP_INPUT_CONFIGURATION_H_ */
