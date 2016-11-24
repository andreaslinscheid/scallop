/*	This file GeneralizedSusceptibility.h is part of scallop.
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
 *  Created on: Nov 11, 2016
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_GW_FLEX_GENERALIZEDSUSCEPTIBILITY_H_
#define SCALLOP_GW_FLEX_GENERALIZEDSUSCEPTIBILITY_H_

#include "scallop/gw_flex/MatsubaraImagTimeFourierTransform.h"
#include "scallop/gw_flex/GreensFunctionOrbital.h"
#include "scallop/gw_flex/MemoryLayout.h"
#include "scallop/gw_flex/InteractionMatrix.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
class GeneralizedSusceptibility : public MatsubaraImagTimeFourierTransform<T>, private MemoryLayout
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;

	using MemoryLayout::get_nOrb;
	using MemoryLayout::get_nChnls;

	GeneralizedSusceptibility();

	void compute_from_gf(GreensFunctionOrbital<T> const & GF,
			size_t channels = 4);

	void initialize_zero(
			size_t nM,
			size_t nO,
			size_t channels,
			std::vector<size_t> grid,
			bool intimeSpace,
			bool inKSpace);

	void spin_RPA_enhancement(
			InteractionMatrix<T> const& interMat);

	T operator() (size_t ik, size_t iw, size_t j, size_t jp, size_t m1,  size_t m2) const;

	T & operator() (size_t ik, size_t iw, size_t j, size_t jp, size_t m1,  size_t m2);

	T operator() (size_t ik, size_t iw, size_t j, size_t jp, size_t l1, size_t l2, size_t l3, size_t l4) const;

	T & operator() (size_t ik, size_t iw, size_t j, size_t jp, size_t l1, size_t l2, size_t l3, size_t l4);

	/**
	 *	Declare that the next compute_from_gf will have to reinitialize.
	 */
	void set_uninitialized();

	auxillary::LinearAlgebraInterface<T> const& get_linAlg_module() const;
private:

	bool isInit_ = false;

	typename auxillary::TemplateTypedefs<T>::scallop_vector bufferQ1_;

	typename auxillary::TemplateTypedefs<T>::scallop_vector bufferQ2_;

	typename auxillary::TemplateTypedefs<T>::scallop_vector bufferSBlock_;

	///We keep a linear algebra module so that it keeps its buffers allocated
	auxillary::LinearAlgebraInterface<T> linAlgModule_;

	void v_matrix_multiplication(
			size_t j, size_t l1, size_t a1, size_t s1, size_t l2,  size_t a2, size_t s2,
			MemoryLayout const& gf_layout,
			size_t &gf_index, T & prefactor) const;
};

} /* namespace gw_flex */
} /* namespace scallop */


#include "scallop/gw_flex/src/GeneralizedSusceptibility.hpp"
#endif /* SCALLOP_GW_FLEX_GENERALIZEDSUSCEPTIBILITY_H_ */
