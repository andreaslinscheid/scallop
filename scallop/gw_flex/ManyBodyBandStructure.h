/*	This file ManyBodyBandStructure.h is part of scallop.
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

#ifndef SCALLOP_GW_FLEX_MANYBODYBANDSTRUCTURE_H_
#define SCALLOP_GW_FLEX_MANYBODYBANDSTRUCTURE_H_

#include "scallop/parallel/IrregularGridDistribution.h"
#include "scallop/input/KPath.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
class ManyBodyBandStructure
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;
	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector V;
	typedef typename auxillary::TemplateTypedefs<bT>::scallop_vector VbT;

	/**
	 * Compute the spectral function of a given data structure on the nearly real axis.
	 *
	 * The method is general so that it can handle both the fermionic Green's function and the susceptiblity.
	 * A template unit '\p blockTransformation' can be used to change the behavior of the method along the
	 * path defined by the irregular grid \p path . Methods 'before' and 'after' allow to perform custom
	 * contribution along the path such as the single particle spectrum and perform post-processing of
	 * the data on the nearly real axis, such a a summation over bands.
	 *
	 * @param path					The distributed among processors irregular points on which the data will be computed.
	 * @param omega					The frequency points, where the data will be computed
	 * @param inInKSpace			If true, the points are in k space.
	 * @param rGrid					The regular grid.
	 * @param regularGridData		The data on the regular grid per processor, ordered according defined by the regular grid.
	 * @param MatsGrid				The Matsubara frequency points.
	 * @param blockD				The number of data elements on a given frequency and space point.
	 * @param blockTransformation	A class that implements
	 * 									0) (optional) a typedef container_type. This can be used to specify in what random acess container
	 * 										intermediate data will processes, e.g. of reduced bit number such as a float.
	 * 										If no container_type is typedefed, the container_type will be of type D.
	 *
	 * 									1) a method 'size_t size_per_block() const' that returns the number data elements for a given
	 * 										frequency and space point.
	 *
	 * 									2) a method 'size_t size_per_block_after() const' that returns the number data elements
	 * 										after the analytic continuation.
	 *
	 * 									3) the method 'before' with the signature
	 *										void before( container_type  const& dataReg, parallel::IrregularGridDistribution<T> const& path,
	 *												V const& MatsGrid) const;
	 *										which allows to add things that are easily obtained at irregualr points and do not have
	 *										to be interpolated such as the single particle spectrum.
	 *										This is applied immediately _before_ the analytic continuation.
	 *
	 * 									4) a method 'after' with the signature
	 * 										void after(container_type const& dataPath, V const& omega, size_t ik,
	 * 											 container_type & mappingData, V & mappingOmega) const;
	 * 										which allows to perform computation on the continued data _after_ the analytic calculations
	 * @param result				The array with the distributed data according to the layout space-pt,frequency,block.
	 */
	template<class B, class D, class DO>
	void compute_spectral_function(
			parallel::IrregularGridDistribution<T> const& path,
			V const& omega,
			bool inInKSpace,
			parallel::GridDistribution<T> const& rGrid,
			D const& regularGridData,
			V const& MatsGrid,
			size_t blockD,
			B const& blockTransformation,
			DO & result);

private:
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/ManyBodyBandStructure.hpp"
#endif /* SCALLOP_GW_FLEX_MANYBODYBANDSTRUCTURE_H_ */
