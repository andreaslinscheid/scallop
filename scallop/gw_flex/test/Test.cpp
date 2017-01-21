/*	This file Test.cpp is part of scallop.
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
 *  Created on: Oct 28, 2016
 *      Author: alinsch
 */

#include "scallop/gw_flex/test/Test.h"
#include "scallop/gw_flex/test/FFTBase_test.h"
#include "scallop/gw_flex/test/MatsubaraImagTimeFourierTransform_test.h"
#include "scallop/gw_flex/test/GreensFunctionOrbital_test.h"
#include "scallop/gw_flex/test/UnitaryWannierKSBands_test.h"
#include "scallop/gw_flex/test/GeneralizedSusceptibility_test.hpp"
#include "scallop/gw_flex/test/PhononGreensFunction_test.hpp"
#include "scallop/gw_flex/test/InteractionMatrix_test.hpp"
#include "scallop/gw_flex/test/SelfEnergy_test.hpp"
#include "scallop/gw_flex/test/KohnShamBandStructure_test.hpp"
#include "scallop/gw_flex/test/ManyBodyBandStructure_test.hpp"
#include "scallop/gw_flex/test/ChemicalPotentialShifting_test.hpp"
#include "scallop/gw_flex/test/DysonEquation_test.hpp"
#include "scallop/gw_flex/test/GapFileReader_test.hpp"
#include <complex>

namespace scallop
{
namespace gw_flex
{
namespace test
{

void Test::run_test()
{
	test_GapFileReader();
	test_selfEnergy();
	test_DysonEquation();
	test_phonon_gf();
	test_InteractionMatrix();
	test_FFTBase();
	test_time_freq_Fourier_transform();
	test_UnitaryWannierKSBands();
	test_GreensFunctionOrbital();
	test_Susceptibility();
	test_KSBandstructure();
	test_MandyBodyBandStructure();
	test_ChemicalPotentialShifting();
}

void Test::test_GapFileReader()
{
	GapFileReader_test g;
	g.test_all();
}

void Test::test_GreensFunctionOrbital()
{
	GreensFunctionOrbital_test< std::complex<double> > gf_test;
	gf_test.test_all();
}

void Test::test_time_freq_Fourier_transform()
{
	MatsubaraImagTimeFourierTransform_test<std::complex<double> > fft_freq_d;
	fft_freq_d.test_free_particle_greensfunction();
}

void Test::test_FFTBase()
{
	typedef std::complex<double> T;
	FFTBase_test<T> fftbaset;
	fftbaset.test_all();
}

void Test::test_UnitaryWannierKSBands()
{
	typedef std::complex<double> T;
	UnitaryWannierKSBands_test<T> unit_tst;
	unit_tst.test_all();
}

void Test::test_Susceptibility()
{
	typedef std::complex<double> T;
	GeneralizedSusceptibility_test<T> stst;
	stst.test_all();
}

void Test::test_phonon_gf()
{
	typedef std::complex<double> T;
	PhononGreensFunction_test<T> phtest;
	phtest.test_all();
}

void Test::test_InteractionMatrix()
{
	typedef std::complex<double> T;
	InteractionMatrix_test<T> Itest;
	Itest.test_all();
}

void Test::test_selfEnergy()
{
	typedef std::complex<double> T;
	SelfEnergy_test<T> setest;
	setest.test_all();
}

void Test::test_KSBandstructure()
{
	typedef std::complex<double> T;
	KohnShamBandStructure_test<T> kstest;
	kstest.test_all();
}

void Test::test_MandyBodyBandStructure()
{
	ManyBodyBandStructure_test< std::complex<double> > mb_test;
	mb_test.test_all();
}

void Test::test_DysonEquation()
{
	DysonEquation_test d_test;
	d_test.test_all();
}

void Test::test_ChemicalPotentialShifting()
{
	ChemicalPotentialShifting_test< std::complex<double> > c;
	c.test_all();
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
