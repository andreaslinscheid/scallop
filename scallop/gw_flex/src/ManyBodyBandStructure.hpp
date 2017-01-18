/*	This file ManyBodyBandStructure.hpp is part of scallop.
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
#include "scallop/auxillary/PadePolynom.h"

namespace scallop
{
namespace gw_flex
{

namespace detail
{

template <class T, class default_Type>
struct defined_container_type
{
private:
    template <typename T1>
    static typename T1::container_type test(int);
    template <typename>
    static default_Type test(...);
public:
    typedef decltype(test<T>(0)) type;
};

} /* namespace detail */

template<typename T>
template<class B, class D, class DO>
void ManyBodyBandStructure<T>::compute_spectral_function(
		parallel::IrregularGridDistribution<T> const& path,
		V const& omega,
		bool inInKSpace,
		parallel::GridDistribution<T> const& rGrid,
		D const& regularGridData,
		V const& MatsGrid,
		size_t blockD,
		B const& blockTransformation,
		DO & result)
{
	typedef typename detail::defined_container_type<B,D>::type container_type;

	size_t nM = MatsGrid.size();
#ifndef NDEBUG
	size_t nGR = inInKSpace ? rGrid.get_num_k_grid() : rGrid.get_num_R_grid();
	assert( nGR*nM*blockD == regularGridData.size() );
	assert( path.get_spaceGrid_proc().get_grid() == rGrid.get_grid() );
#endif

	container_type dataForThisPts;
	path.proc_sync_data(inInKSpace, regularGridData,  dataForThisPts, nM*blockD );

	container_type dataAlongPath;
	path.linear_interpolate_data(dataForThisPts,dataAlongPath,nM*blockD);

	blockTransformation.before( dataAlongPath, path, MatsGrid );

	typename auxillary::TemplateTypedefs< std::complex<double> >::scallop_vector thisFreqBlock(nM);
	result = typename std::remove_reference<DO>::type( blockD*omega.size()* path.get_n_pts() );
	for ( size_t ikPath = 0; ikPath < path.get_n_pts(); ++ikPath )
	{
		assert(dataAlongPath.size() >= path.get_n_pts()*blockD*nM);

		for ( size_t iB = 0 ; iB < blockD; ++iB)
		{
			for ( size_t iM = 0 ; iM < nM; ++iM)
				thisFreqBlock[iM] = dataAlongPath[(ikPath*nM+iM)*blockD+iB];

			auxillary::PadePolynom fOnRealAxis( MatsGrid.data(), thisFreqBlock.data(), nM );
			for (size_t io = 0 ; io< omega.size(); ++io )
				result[(ikPath*omega.size()+io)*blockD+iB] = fOnRealAxis( omega[io] );
		}
	}
	blockTransformation.after( result, path, omega);
}

} /* namespace gw_flex */
} /* namespace scallop */
