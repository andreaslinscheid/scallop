/*	This file SusceptibilityPlotter.h is part of scallop.
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
 *  Created on: Feb 7, 2017
 *      Author: A. Linscheid
 */

#ifndef SCALLOP_OUTPUT_SUSCEPTIBILITYPLOTTER_H_
#define SCALLOP_OUTPUT_SUSCEPTIBILITYPLOTTER_H_

#include "scallop/output/LowestMatsFreqOnSpaceGridPlotter.h"
#include "scallop/auxillary/TemplateTypedefs.h"
#include "scallop/output/Projector3DTo2D.h"
#include "scallop/parallel/GridDistribution.h"

namespace scallop
{
namespace output
{

class SusceptibilityPlotter : public LowestMatsFreqOnSpaceGridPlotter
{
public:
	/**
	 * Calling this means constructor means do not plot
	 */
	SusceptibilityPlotter();

	void configure_spin(bool plotSpinSus);

	void configure_charge(bool plotChargeSus);

	void set_plot_physical_sus();

	void set_plot_largest_eval();

	SusceptibilityPlotter( std::vector<size_t> kgrid,
			size_t NOrbitals, size_t bufferDim, size_t channels);

	void set_buffer(std::vector<size_t> kgrid,
			size_t NOrbitals, size_t bufferDim, size_t channels);

	template<class D>
	void append_data_to_buffer(D const& data);

	void plot_static_k();

	bool do_plot_static() const;

	bool do_plot_spin() const;

	bool do_plot_charge() const;

	typedef std::vector<size_t> (*Projector) (std::vector<size_t> const & kvector);

	void define_projection_plane_3D(Projector3DTo2D p);
private:

	bool plotSpin_ = false;

	bool plotCharge_ = false;

	bool plotPhysicalSus_ = false;

	bool plotLargestEVal_ = true;

	parallel::GridDistribution<float> kgrid_;

	size_t nO_ = 0;

	size_t nC_ = 0;

	auxillary::TemplateTypedefs< std::complex<float> >::scallop_vector buffer_;

	size_t pos_ = 0;

	Projector3DTo2D p_;
};

} /* namespace output */
} /* namespace scallop */

#include "scallop/output/src/SusceptibilityPlotter.hpp"
#endif /* SCALLOP_OUTPUT_SUSCEPTIBILITYPLOTTER_H_ */
