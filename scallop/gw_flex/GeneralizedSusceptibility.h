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
#include "scallop/output/SusceptibilityPlotter.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
class GeneralizedSusceptibility : public MatsubaraImagTimeFourierTransform<T>, private MemoryLayout
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;

	/**
	 * Control the adiabatic switching on of the interaction parameter in case of an instability.
	 *
	 * Passing an instance of this object to RPA_enhancement
	 * allows to modify the attempts made to circumvent an instability
	 * by reducing the interaction homogeneously up to maxScaling.
	 * The code will try scaleResolutionSteps to reduce the scaling each time
	 * halving the interval from scaling to maxScaling starting with scaling=1.0.
	 * Given the instability is avoided if the interaction is scaled down by maxScaling,
	 * the code will try maxScaleSteps with the above procedure before giving up and
	 * collecting the smallest eigenvalue and vector.
	 */
	typedef class AdiabaticUpscale_
	{
	public:

		void set(bT maxScaling, size_t maxScaleResolutionSteps, size_t maxScaleSteps);

		bool is_soft() const;

		bool steps_ok() const;

		void reset();

		void get_soft_vector( std::vector<bT> & softKVector, std::vector<T> & softOrbitalVector) const;

		void track_instability(  std::vector<bT> const& softKVector, std::vector<T> const& softOrbitalVector );

		bT get_scaling() const;

		//Resets the scaling, but if currently, the scaling was not 1, then we increase the nscale counter,
		//since this means a stable susceptiblity was only possible when reducing the interaction.
		void set_stable();
	private:

		bT maxScaling_ = bT(1);

		bT scaling_ = bT(1);

		size_t maxScaleResolutionSteps_ = 1;

		size_t nscaleResolutionSteps_ = 0;

		size_t nscale_ = 0;

		size_t maxNScale_ = 0;

		std::vector<bT> softKVector_;

		std::vector<T> softOrbitalVector_;

		bT softEV_ = bT(0);

	} AdiabaticUpscale;

	using MemoryLayout::get_nOrb;
	using MemoryLayout::get_nChnls;

	GeneralizedSusceptibility();

	void compute_from_gf(
			GreensFunctionOrbital<T> const & GF,
			size_t channels,
			size_t channel_offset);

	void plot_static(
			output::SusceptibilityPlotter & plotter) const;

	void initialize_zero(
			size_t nM,
			size_t nO,
			size_t channels,
			std::vector<size_t> grid,
			bool intimeSpace,
			bool inKSpace);

	auxillary::LinearAlgebraInterface<T> const& get_linAlg_module() const;

	void RPA_enhancement(
			InteractionMatrix<T> const& interMat,
			AdiabaticUpscale & a,
			output::SusceptibilityPlotter & plotter );

	void RPA_enhancement(
			InteractionMatrix<T> const& interMat,
			AdiabaticUpscale & a );

	void RPA_enhancement(
			InteractionMatrix<T> const& interMat,
			output::SusceptibilityPlotter & plotter );

	void RPA_enhancement(
			InteractionMatrix<T> const& interMat);

	bool has_charge_part() const;

	size_t channels_offset() const;

protected:
	T operator() (size_t ik, size_t iw, size_t j, size_t jp, size_t m1,  size_t m2) const;

	T & operator() (size_t ik, size_t iw, size_t j, size_t jp, size_t m1,  size_t m2);

	T operator() (size_t ik, size_t iw, size_t j, size_t jp, size_t l1, size_t l2, size_t l3, size_t l4) const;

	T & operator() (size_t ik, size_t iw, size_t j, size_t jp, size_t l1, size_t l2, size_t l3, size_t l4);

private:

	typename auxillary::TemplateTypedefs<T>::scallop_vector bufferQ1_;

	typename auxillary::TemplateTypedefs<T>::scallop_vector bufferQ2_;

	typename auxillary::TemplateTypedefs<T>::scallop_vector bufferSBlock_;

	///We keep a linear algebra module so that it keeps its buffers allocated
	auxillary::LinearAlgebraInterface<T> linAlgModule_;

	///Determines if we compute the pure charge and spin+charge (0) or the pure spin suceptibility (1)
	size_t channelsOffset_ = 0;

	void v_matrix_multiplication(
			size_t j, size_t l1, size_t a1, size_t s1, size_t l2,  size_t a2, size_t s2,
			MemoryLayout const& gf_layout,
			size_t &gf_index, T & prefactor) const;
};

} /* namespace gw_flex */
} /* namespace scallop */


#include "scallop/gw_flex/src/GeneralizedSusceptibility.hpp"
#endif /* SCALLOP_GW_FLEX_GENERALIZEDSUSCEPTIBILITY_H_ */
