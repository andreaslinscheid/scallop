/*	This file DysonEquation_test.hpp is part of scallop.
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
 *  Created on: Dec 16, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_DYSONEQUATION_TEST_HPP_
#define SCALLOP_GW_FLEX_DYSONEQUATION_TEST_HPP_

#include "scallop/gw_flex/DysonEquation.h"
#include "scallop/output/TerminalOut.h"
#include "scallop/auxillary/BasicFunctions.h"
#include "scallop/gw_flex/KohnShamGreensFunctionOrbital.h"
#include "scallop/gw_flex/KohnShamBandStructure.h"
#include "scallop/gw_flex/test/KohnShamBandStructure_test.hpp"
#include "scallop/gw_flex/GreensFunctionOrbital.h"

namespace scallop
{
namespace gw_flex
{
namespace test
{

class DysonEquation_test
{
public:

	void test_all();

private:

	void test_zero_se();
};

void DysonEquation_test::test_all()
{
	this->test_zero_se();
}

void DysonEquation_test::test_zero_se()
{
	typedef std::complex<double> T;

	const double temp = 10;
	const size_t nM = 32;
	const double beta = auxillary::BasicFunctions::inverse_temperature( temp );

	//Construct a KS Greens function from the solution of the Dyson equation with zero SE
	output::TerminalOut msg;
	msg << "Testing the Dyson equation module by construction the KS GF.";

	msg << "Computing the reference KS GF:";
	std::vector<size_t> grid = {16, 16};
	KohnShamBandStructure<T> ksbnd;
	const std::string filename = "/tmp/2bndcos.dat";
	KohnShamBandStructure_test<T> ksBndTest;
	ksBndTest.write_test_input_file( filename );
	ksbnd.initialize_from_file( grid , filename );
	KohnShamGreensFunctionOrbital<T> ksgfRef;
	ksgfRef.set_from_KS_bandstructure(false,nM, beta, ksbnd);

	msg << "Solving the Dyson equation:";
	DysonEquation d;
	GreensFunctionOrbital<T> g = ksgfRef;
	SelfEnergy<T> se;
	typename auxillary::TemplateTypedefs<T>::scallop_vector data;
	se.initialize( nM, grid, g.get_nOrb(), false, true, data, ksbnd.get_chem_pot() );

	auto m = auxillary::BasicFunctions::matzubara_frequency_array(nM,beta);
	d.solve_by_inversion(g,m ,ksbnd, se);

	for ( size_t ik = 0 ; ik < g.get_spaceGrid_proc().get_num_k_grid(); ++ik)
	{
		for ( size_t iw = 0 ; iw < nM; ++iw)
		{
			auto ptr_g = g.read_data_ptr_block(ik,iw);
			auto ptr_gks = ksgfRef.read_data_ptr_block(ik,iw);
			for ( size_t ib = 0 ; ib < g.get_data_block_size(); ++ib)
			{
				const double maxDiff = 1e-6;
				if ( std::abs(ptr_g[ib]-ptr_gks[ib]) > maxDiff )
					error_handling::Error("Test failed at (ik,iw,ib)=("
							+std::to_string( ik )+","+std::to_string( iw )+","+std::to_string( ib )+
							"), reference KS GF ("
							+std::to_string( std::real(ptr_gks[ib]) )+","+std::to_string( std::imag(ptr_gks[ib]) )+
							") and Dyson equation solution ("
							+std::to_string( std::real(ptr_g[ib]) )+","+std::to_string( std::imag(ptr_g[ib]) )+
							") differ by ("
							+std::to_string( std::real(ptr_g[ib]-ptr_gks[ib]) )+","+std::to_string( std::imag(ptr_g[ib]-ptr_gks[ib]) )+
							"), in abs more than "+std::to_string(maxDiff));
			}
		}
	}
	msg <<"done.";
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */

#endif /* SCALLOP_GW_FLEX_DYSONEQUATION_TEST_HPP_ */
