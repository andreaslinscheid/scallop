/*	This file SusceptibilityPlotter.cpp is part of scallop.
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
 *  Created on: Feb 8, 2017
 *      Author: A. Linscheid
 */

#include "scallop/output/SusceptibilityPlotter.h"
#include "scallop/auxillary/LinearAlgebraInterface.h"
#include "scallop/gw_flex/GeneralizedSusceptibility.h"
#include "scallop/parallel/IrregularGridDistribution.h"

namespace scallop
{
namespace output
{

SusceptibilityPlotter::SusceptibilityPlotter()
{

}

SusceptibilityPlotter::SusceptibilityPlotter(
		std::vector<size_t> kgrid,
		size_t NOrbitals,
		size_t bufferDim,
		size_t channels)
{
	this->set_buffer(kgrid, NOrbitals, bufferDim, channels);
}


void SusceptibilityPlotter::set_plot_physical_sus()
{
	plotLargestEVal_ = false;
	plotPhysicalSus_ = true;
}

void SusceptibilityPlotter::set_plot_largest_eval()
{
	plotLargestEVal_ = true;
	plotPhysicalSus_ = false;
}

void SusceptibilityPlotter::configure_spin(bool plotSpinSus)
{
	plotSpin_ = plotSpinSus;
}

void SusceptibilityPlotter::configure_charge(bool plotChargeSus)
{
	plotCharge_ = plotChargeSus;
}

void SusceptibilityPlotter::set_buffer(std::vector<size_t> kgrid,
		size_t NOrbitals, size_t bufferDim, size_t channels)
{
	kgrid_.distribute_grid( kgrid );
	nO_ = NOrbitals;
	nC_ = channels;
	buffer_ = auxillary::TemplateTypedefs< std::complex<float> >::scallop_vector(bufferDim, std::complex<float>(0) );
	pos_ = 0;
}

bool SusceptibilityPlotter::do_plot_static() const
{
	return kgrid_.get_grid().size() > 0;
}

bool SusceptibilityPlotter::do_plot_spin() const
{
	return plotSpin_;
}

bool SusceptibilityPlotter::do_plot_charge() const
{
	return plotCharge_;
}

void SusceptibilityPlotter::plot_static_k(  )
{
	assert( kgrid_.get_grid().size() > 0 );

	typedef std::complex<float> T;
	size_t nK = 1;
	for ( auto kx : kgrid_.get_grid() )
		nK *= kx;

	size_t B = nO_*nO_*nC_;
	size_t B2 = B*B;
	assert( B2*nK == buffer_.size() );

	std::vector< float > reducedSust( nK );
	if ( plotLargestEVal_)
	{
		auxillary::LinearAlgebraInterface<T> linalg;
		std::vector<float> ev(B);
		for ( size_t ik = 0 ; ik < nK; ++ik )
		{
			linalg.hermitian_eigensystem(true, false, &( buffer_[ik*B2] ),ev.size(),ev.data() );
			reducedSust[ik] = *std::min_element(ev.begin(),ev.end());
		}
	}
	else if ( plotPhysicalSus_ )
	{
		gw_flex::MemoryLayout mem;
		mem.initialize_layout_4pt_scalar_obj( nO_, nC_ );
		for ( size_t ik = 0 ; ik < nK; ++ik )
		{
			decltype(buffer_)::value_type tmp = 0;
			for ( size_t l1=0;  l1 < nO_; ++l1 )
				for ( size_t l2=0;  l2 < nO_; ++l2 )
				{
					tmp += buffer_[ik*B2+mem.memory_layout_4pt_scalar_obj(0,0,l1,l1,l2,l2)];
				}
			reducedSust[ik] = std::real(tmp);
			//make sure the static susceptibility is real
			assert( std::abs(std::imag(tmp)) < 0.000001 );
		}
	}
	else
	{
		error_handling::Error("No plot property defined");
	}

	//Check if the final plot is 2D
	if ( ((kgrid_.get_grid().size() == 2) or ((kgrid_.get_grid().size() == 3) && p_.project_active() ))
			and (this->get_plot_program() == PlotProgram::defaultPlotProg) )
	{
		this->set_active_plotprogram( PlotProgram::Gnuplot );
	}

	if ( kgrid_.get_grid().size() == 2 )
	{
		std::vector<float> kx(kgrid_.get_grid()[0]);
		for ( size_t ix = 0 ; ix < kx.size(); ++ix )
			kx[ix] = float(ix)/float(kx.size());
		std::vector<float> ky(kgrid_.get_grid()[1]);
		for ( size_t iy = 0 ; iy < ky.size(); ++iy )
			ky[iy] = float(iy)/float(ky.size());

		this->plot_2D_grid_data(reducedSust,ky,kx);
	}

	//3D data but projection plane
	if ( (kgrid_.get_grid().size() == 3)  && p_.project_active() )
	{
		auto kgrid2D = p_.get_2D_grid_dimensions();

		//Obtain the grid points for gnuplot
		std::vector<float> kx(kgrid2D[0]);
		for ( size_t ix = 0 ; ix < kx.size(); ++ix )
			kx[ix] = float(ix)/float(kx.size());
		std::vector<float> ky(kgrid2D[1]);
		for ( size_t iy = 0 ; iy < ky.size(); ++iy )
			ky[iy] = float(iy)/float(ky.size());

		//define the grid in the order as the parallel grid distribution
		//perform the mapping and associate the data
		parallel::GridDistribution<float> PlotGrid;
		PlotGrid.distribute_grid( kgrid2D );
		auxillary::TemplateTypedefs<float>::scallop_vector allPoints(kx.size()*ky.size()*3);
		for (size_t ikplot = 0 ; ikplot < PlotGrid.get_num_grid(); ++ikplot)
		{
			auto kVectorPlot = PlotGrid.k_conseq_to_xyz(ikplot);
			auto mappedk = p_(kVectorPlot);
			assert( mappedk.size() == 3 );
			for ( size_t i = 0 ; i < 3 ; ++i)
				allPoints[ikplot*3+i] = float(mappedk[i])/float(kgrid_.get_grid()[i]);
		}
		parallel::IrregularGridDistribution<float> twoDgrid;
		twoDgrid.distribute_pts( true, allPoints, kgrid_ );

		//scatter the data into the plot grid
		std::vector< float > dataProcPlot;
		twoDgrid.proc_sync_data( true, reducedSust, dataProcPlot, 1);

		std::vector< float > dataOnPlotGrid;
		twoDgrid.linear_interpolate_data(dataProcPlot,dataOnPlotGrid,1);

		this->plot_2D_grid_data(dataOnPlotGrid,ky,kx);
	}
}

void SusceptibilityPlotter::define_projection_plane_3D(Projector3DTo2D p)
{
	p_ = p;
}

} /* namespace output */
} /* namespace scallop */
