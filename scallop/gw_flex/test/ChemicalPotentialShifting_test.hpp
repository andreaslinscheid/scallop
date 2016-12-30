/*	This file ChemicalPotentialShifting_test.hpp is part of scallop.
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
 *  Created on: Dec 17, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_TEST_CHEMICALPOTENTIALSHIFTING_TEST_HPP_
#define SCALLOP_GW_FLEX_TEST_CHEMICALPOTENTIALSHIFTING_TEST_HPP_

#include "scallop/gw_flex/ChemicalPotentialShifting.h"
#include "scallop/gw_flex/KohnShamGreensFunctionOrbital.h"
#include "scallop/gw_flex/GreensFunctionOrbital.h"
#include "scallop/gw_flex/KohnShamBandStructure.h"
#include "scallop/auxillary/TemplateTypedefs.h"
#include "scallop/auxillary/TypeMapComplex.h"
#include "scallop/auxillary/BasicFunctions.h"
#include "scallop/gw_flex/test/KohnShamBandStructure_test.hpp"

namespace scallop
{
namespace gw_flex
{
namespace test
{

template<typename T>
class ChemicalPotentialShifting_test
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;

	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector V;

	void test_all();

private:

	void shift_cosine();

};

template<typename T>
void ChemicalPotentialShifting_test<T>::test_all()
{
	this->shift_cosine();
}

template<typename T>
void ChemicalPotentialShifting_test<T>::shift_cosine()
{
	output::TerminalOut msg;

	msg << "Testing the chemical potential adjustment for a two band cosine model.";

	std::pair<bool,bT> result;
	const bT temp = 10;
	const bT Ne = 3;
	auto beta = auxillary::BasicFunctions::inverse_temperature( temp );
	std::vector<size_t> grid = {32,32};
	KohnShamBandStructure_test<T> bndtest;
	const std::string filename = "/tmp/2bndcos.dat";
	bndtest.write_test_input_file(filename);

	ChemicalPotentialShifting<T> chmpot;

	KohnShamBandStructure<T> elstr;
	elstr.initialize_from_file( grid, filename );
	elstr.adjust_filling( Ne, beta );
	auto NeC = elstr.compute_N_electrons( beta );
	if ( std::abs( NeC - Ne) > 0.01 )
		error_handling::Error(std::string("Error adjusting number of electrons, not ")
				+std::to_string(Ne)+" but "+std::to_string(NeC));

	size_t nM = 64;

	KohnShamGreensFunctionOrbital<T> ksGF;
	ksGF.set_from_KS_bandstructure( false, nM, beta, elstr );
	ksGF.perform_space_fft();
	GreensFunctionOrbital<T> g = ksGF;

	chmpot.determine_new_chemPot( Ne, beta, elstr, ksGF, g, result);
	if ( ! result.first )
	{
		error_handling::Error("Test failed, chemical potential shifting did not succeed");
	}
	elstr.set_chem_pot( result.second );
	ksGF.set_chem_pot( result.second );
	g.set_chem_pot( result.second );

	msg << "done.";
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */

#endif /* SCALLOP_GW_FLEX_TEST_CHEMICALPOTENTIALSHIFTING_TEST_HPP_ */
