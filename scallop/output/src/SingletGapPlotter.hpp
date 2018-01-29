/*	This file SingletGapPlotter.hpp is part of scallop.
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
 *  Created on: Feb 6, 2017
 *      Author: A. Linscheid
 */

#include "scallop/output/SingletGapPlotter.h"
#include "scallop/gw_flex/MemoryLayout.h"

namespace scallop
{
namespace output
{

template<typename T>
void SingletGapPlotter::plot(PlotDataProcessing<T> data)
{
	//get the pair data
	gw_flex::SelfEnergy<T> const & se = data.get_self_energy();

	assert( se.is_in_k_space() && ( not se.is_in_time_space() ) );
	size_t nK = se.get_spaceGrid_proc().get_num_k_grid();
	size_t nO = se.get_nOrb();

	gw_flex::MemoryLayout mem;
	mem.initialize_layout_2pt_obj(nO);

	std::vector< std::complex<float> > pairData( nO*nO*nK );
	for ( size_t ik = 0 ; ik < nK; ++ik )
	{
		auto d_ptr = se.read_phs_grid_ptr_block(ik,0);
		for ( size_t l1 = 0; l1 < nO; ++l1)
			for ( size_t l2 = 0; l2 < nO; ++l2)
			{
				pairData[(ik*nO+l1)*nO+l2] =
					d_ptr[mem.memory_layout_2pt_obj(l1,0,0,l2,1,1)]
					-d_ptr[mem.memory_layout_2pt_obj(l1,0,1,l2,1,0)]
					-d_ptr[mem.memory_layout_2pt_obj(l1,1,0,l2,0,1)]
					+d_ptr[mem.memory_layout_2pt_obj(l1,1,1,l2,0,0)];
			}

	}

	if ( se.get_spaceGrid_proc().get_dim() == 2 )
	{
		if ( this->get_plot_program() == PlotProgram::defaultPlotProg )
		{
			this->set_active_plotprogram( PlotProgram::Gnuplot );
		}

		std::vector<float> kx(se.get_spaceGrid_proc().get_grid()[0]);
		for ( size_t ix = 0 ; ix < kx.size(); ++ix )
			kx[ix] = float(ix)/float(kx.size());
		std::vector<float> ky(se.get_spaceGrid_proc().get_grid()[1]);
		for ( size_t iy = 0 ; iy < ky.size(); ++iy )
			ky[iy] = float(iy)/float(ky.size());

		//todo Here we do a dirty hack ...
		std::vector< float > pairPrj( nK );

		for ( size_t ik = 0 ; ik < nK; ++ik )
			pairPrj[ik] = std::real(pairData[(ik*nO+0)*nO+0]);
		this->plot_2D_grid_data(pairPrj,ky,kx);
	}
}

} /* namespace output */
} /* namespace scallop */
