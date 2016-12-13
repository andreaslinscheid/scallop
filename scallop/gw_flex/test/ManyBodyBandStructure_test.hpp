/*	This file ManyBodyBandStructure_test.hpp is part of scallop.
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
 *  Created on: Nov 29, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/ManyBodyBandStructure.h"
#include "scallop/output/TerminalOut.h"
#include "scallop/gw_flex/test/GreensFunctionOrbital_test.h"
#include "scallop/input/KPath.h"
#include "scallop/parallel/IrregularGridDistribution.h"
#include "scallop/auxillary/BasicFunctions.h"
#include "scallop/gw_flex/ProjectOutBandStructure.h"

namespace scallop
{
namespace gw_flex
{
namespace test
{

template<typename T>
class ManyBodyBandStructure_test
{
public:

	void test_all();
private:

	void test_construction();
};

template<typename T>
void ManyBodyBandStructure_test<T>::test_all()
{
	test_construction();
}

template<typename T>
void ManyBodyBandStructure_test<T>::test_construction()
{
	typedef typename ManyBodyBandStructure<T>::V V;

	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	typedef typename auxillary::TypeMapComplex<T>::type bT;
	output::TerminalOut msg;

	msg << "Testing the ManyBodyBandStructure unit on a cosine band model.";

	std::vector<size_t> spaceGrid = { 64, 64 };
	const bT bandwidth = 20; // meV
	const bT mu = 0; // Perfect nesting
	const bT temperature = 0.0376875	; // Temperature hand tuned to yield -0.1 to good accuracy at Q=pi,pi
	const bT kb = 0.86173324 ; // meV / K
	const bT beta = 1.0 / (kb * temperature);
	const size_t nM = 128;

	auto cos_bnd = [&] (std::vector<size_t> const& tuple,
			size_t l1, size_t a1, size_t s1,
			size_t l2, size_t a2, size_t s2){
		return ( a1 == a2 ? a1 == 0 ? 1.0 : -1.0 : 0.0 ) *  ( s1 == s2 ? 1.0 : 0.0 )
				*( 0.5*bandwidth*(std::cos( (2*M_PI*tuple[0])/spaceGrid[0] )
							+std::cos( (2*M_PI*tuple[1])/spaceGrid[1] ))-mu);
	};

	GreensFunctionOrbital_test<T> gfTestUnit;
	auto gf = gfTestUnit.construct_gf_bnd(
			temperature, spaceGrid ,
			nM, 1,
			cos_bnd);


	std::string const & kpathFileName = "/tmp/test_path.dat";
	if ( mpi.ioproc() )
	{
		std::ofstream file(kpathFileName.c_str());
		if ( ! file.good() )
			error_handling::Error(std::string("Unable to open file for writing: ")+kpathFileName);

		file << "X        0.5000  0.0000 100     $\\Gamma$  0.0000  0.0000\n"
			 << "$\\Gamma$  0.0000  0.0000 200     M        0.5000  0.5000\n"
			 << "M        0.5000  0.5000 99      X        0.5000  0.0000\n";
		file.close();
	}

	input::KPath<bT> path;
	path.read_kpath_file(kpathFileName);

	parallel::IrregularGridDistribution<T> irrPath;
	irrPath.distribute_pts(
			/*in k space=*/ true,
			path.get_k_path(),
			gf.get_spaceGrid_proc());

	bT const feqMin = -100;
	bT const feqMax =  100;
	size_t const nOm = 200;
	V omega(nOm,T(0,0.1));
	for (size_t io = 0 ;io  <nOm ; ++io )
		omega[io] += feqMin + (bT(io)+0.5)/bT(nOm)*(feqMax-feqMin);

	V MatsFreqs(nM/2);
	for ( int n = static_cast<int>(nM/2) ; n < static_cast<int>(nM) ; ++n)
	{
		int frequencyIndex =
			(n < static_cast<int>(nM)/2 ? n : n-static_cast<int>(nM) );
		MatsFreqs[n-static_cast<int>(nM/2)] =  T(0,M_PI / beta * ( 2*frequencyIndex+1 ) );
	}

	struct BlockTransformation
	{
		size_t size_per_block() const {return 1;};

		size_t size_per_block_after() const {return 1;};


		void before( V & dataPath, parallel::IrregularGridDistribution<T> const& path, V const& MatsFreq) const
		{
			//we only perform a linear interpolation of everything in this example.
		}

		void after( V & dataPath, parallel::IrregularGridDistribution<T> const& path, V const& omega) const
		{
		};
	};

	V dataOnRegularGrid( gf.get_spaceGrid_proc().get_num_k_grid()*nM/2);
	for ( size_t ik = 0 ; ik < gf.get_spaceGrid_proc().get_num_k_grid(); ++ik )
		for ( size_t iw = 0 ; iw < nM/2; ++iw )
			dataOnRegularGrid[ik*nM/2+iw] = gf.read_phs_grid_ptr_block(ik,iw+nM/2)[0];

	BlockTransformation btrans;

	V result;
	ManyBodyBandStructure<T> mb;
	mb.compute_spectral_function( irrPath,
			omega,
			true,
			gf.get_spaceGrid_proc(),
			dataOnRegularGrid,
			MatsFreqs,
			1,
			btrans,
			result);

	path.band_structure_gnuplot( result, omega, irrPath, "/tmp/test_bnd","$\\varepsilon_{\\boldsymbol{k},n}$ [meV]");
	msg << "done. Data and plotscript at '/tmp/test_bnd.dat' and '/tmp/test_bnd.gp', respectively.";

	auxillary::BasicFunctions bf;
	auto scallopPath = bf.get_scallop_path();

	msg << "\nTesting the ManyBodyBandStructure unit on 2 band consine model.";
	std::string const & modelFileName =scallopPath+"examples/twoBndCosine/model_test.dat";
	if ( mpi.ioproc() )
	{
		std::ofstream file(modelFileName.c_str());
		if ( ! file.good() )
			error_handling::Error(std::string("Unable to open file for writing: ")+modelFileName);

		file << "model TwoBandCosine\n"
				"t1=100\n"
				"t2=150\n"
				"Eh=-5.5\n"
				"Ee=-40\n";
		file.close();
	}

	path.read_kpath_file(scallopPath+"examples/lifeas/input_data/kpath.dat");

	std::vector<size_t> grid = {8,8,4};
	KohnShamBandStructure<T> ksBnd;
	ksBnd.initialize_from_file(grid,modelFileName);

	irrPath.distribute_pts(
			/*in k space=*/ true,
			path.get_k_path(),
			ksBnd.get_unitary().get_spaceGrid_proc());

	V data;
	SelfEnergy<T> emptySE;
	emptySE.initialize(
			nM,
			grid,
			ksBnd.get_nOrb(),
			/* initialInTimeDomain*/ false,
			/* initialInReciprocalDomain */true,
			data);

	MatsFreqs = V(nM);
	for ( int n = 0 ; n < static_cast<int>(nM) ; ++n)
	{
		int frequencyIndex =
			(n < static_cast<int>(nM)/2 ? n : n-static_cast<int>(nM) );
		MatsFreqs[n] =  T(0,M_PI / beta * ( 2*frequencyIndex+1 ) );
	}

	ProjectOutBandStructure<T> bandProj( ksBnd );

	V mappedMatsFreq;
	typename auxillary::TemplateTypedefs<std::complex<float> >::scallop_vector mappedData;
	bandProj.project_out_KS_bands_half_freq(
			emptySE,
			MatsFreqs,
			mappedMatsFreq, mappedData);

	bT const feqMinLi = -150;
	bT const feqMaxLi =  150;
	for (size_t io = 0 ;io  <nOm ; ++io )
		omega[io] = T(feqMinLi + (bT(io)+0.5)/bT(nOm)*(feqMaxLi-feqMinLi),0.01);

	typename ProjectOutBandStructure<T>::container_type result2Cos;
	mb.compute_spectral_function( irrPath,
				omega,
				true,
				emptySE.get_spaceGrid_proc(),
				mappedData,
				mappedMatsFreq,
				emptySE.get_data_block_size(),
				bandProj,
				result2Cos);

	std::string cosinebands = scallopPath+"examples/twoBndCosine/output_data/bnd";
	path.band_structure_gnuplot( result2Cos, omega, irrPath, cosinebands,"$\\varepsilon_{\\boldsymbol{k},n}$ [meV]");
	msg << "done. Data and plotscript at "<<cosinebands;

	msg << "\nTesting the ManyBodyBandStructure unit on 10 band LiFeAs fit.";
	size_t nMLiFeAs = 128;
	grid = {4,4,2};
	ksBnd.initialize_from_file(grid,scallopPath+"examples/lifeas/input_data/LiFeAs_hr.dat");
	ksBnd.set_chem_pot( 9570.0 );

	ProjectOutBandStructure<T> bandProjLiFeAs( ksBnd );

	emptySE.initialize(
			nMLiFeAs,
			grid,
			ksBnd.get_nOrb(),
			/* initialInTimeDomain*/ false,
			/* initialInReciprocalDomain */true,
			data);

	const bT betaLiFeAs = 1.0 / (kb * 10.0);
	MatsFreqs = V(nMLiFeAs);
	for ( int n = 0 ; n < static_cast<int>(nMLiFeAs) ; ++n)
	{
		int frequencyIndex =
			(n < static_cast<int>(nMLiFeAs)/2 ? n : n-static_cast<int>(nMLiFeAs) );
		MatsFreqs[n] =  T(0,M_PI / betaLiFeAs * ( 2*frequencyIndex+1 ) );
	}

	mappedData = typename auxillary::TemplateTypedefs<std::complex<float> >::scallop_vector();
	bandProjLiFeAs.project_out_KS_bands_half_freq(
			emptySE,
			MatsFreqs,
			mappedMatsFreq, mappedData);

	bT const feqMinLiFeAs = -4000;
	bT const feqMaxLiFeAs = 3000;
	for (size_t io = 0 ;io  <nOm ; ++io )
		omega[io] = T(feqMinLiFeAs + (bT(io)+0.5)/bT(nOm)*(feqMaxLiFeAs-feqMinLiFeAs),0.5);

	irrPath.distribute_pts(
			/*in k space=*/ true,
			path.get_k_path(),
			emptySE.get_spaceGrid_proc());

	typename ProjectOutBandStructure<T>::container_type resultLiFeAs;
	mb.compute_spectral_function( irrPath,
				omega,
				true,
				emptySE.get_spaceGrid_proc(),
				mappedData,
				mappedMatsFreq,
				emptySE.get_data_block_size(),
				bandProjLiFeAs,
				resultLiFeAs);

	std::string lifeasbands = scallopPath+"examples/lifeas/input_data/LiFeAs_bnd";
	path.band_structure_gnuplot( resultLiFeAs, omega, irrPath, lifeasbands,"$\\varepsilon_{\\boldsymbol{k},n}$ [meV]");
	msg << "done. Data and plotscript at "<<lifeasbands;
}

} /* namespace test */
} /* namespace gw_flex */
} /* namespace scallop */
