/*	This file LowestMatsFreqOnSpaceGridPlotter.cpp is part of scallop.
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
 *  Created on: Jan 21, 2017
 *      Author: A. Linscheid
 */

#include "scallop/output/LowestMatsFreqOnSpaceGridPlotter.h"
#include "scallop/output/GnuplotBinaryData2DFile.h"
#include <regex>
#include <fstream>

namespace scallop
{
namespace output
{

void LowestMatsFreqOnSpaceGridPlotter::plot_2D_grid_data(
		std::vector<float> data,
		std::vector<float> gridInY,
		std::vector<float> gridInX) const
{
	assert( this->get_plot_program() == FormatDescriptor::PlotProgram::Gnuplot);

	auto dataFile = this->get_filename()+".dat";
	GnuplotBinaryData2DFile bf;
	bf.write_gnuplot_binary_2d_data_file(
			dataFile,
			data,
			gridInY,
			gridInX);

	//The io-proc writes the plot-script
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	auto mm = std::minmax_element(data.begin(), data.end() );
	auto minAll = *mm.first;
	auto maxAll = *mm.second;
	mpi.min( minAll );
	mpi.max( maxAll );

	if ( mpi.ioproc() )
	{
		std::string xTics;
		for ( auto x : gridInX)
			xTics += std::string("\" \"")+std::to_string(float(x))+",";
		xTics.pop_back();

		std::string yTics;
		for ( auto y : gridInY)
			yTics += std::string("\" \"")+std::to_string(float(y))+",";
		yTics.pop_back();

		auto qtty_str = this->get_quantity_label();
		auto quantityLabel = std::regex_replace(qtty_str, std::regex("\\\\"), "\\\\");

		auto ceil_to_digits = [] (double value, int digits)
		{
		    if (value == 0.0)
		        return 0.0;

		    double factor = std::pow(10.0, digits - std::ceil(std::log10(std::fabs(value))));
		    return std::ceil(value * factor) / factor;
		};

		auto inc = ceil_to_digits((maxAll-minAll)/5.0,1);
		std::string colormap =
				"set colorbox horizontal\n"
				"set colorbox user origin 0.25,0.1 size 0.5,0.05\n"
				"set cblabel \""+quantityLabel+"\"\n"
				"set cbtics "+std::to_string(ceil_to_digits(minAll,1))+","+std::to_string(inc)+"\n";

		if ( minAll*maxAll < 0 )
		{
			colormap +=
			"set palette model RGB functions r(gray), g(gray), b(gray)\n";
			"minVal="+std::to_string(minAll)+"\n"
			"maxVal="+std::to_string(maxAll)+"\n"
			"belowZeroInterval=-minVal/abs(maxVal-minVal)\n"
			"aboveZeroInterval=maxVal/abs(maxVal-minVal)\n\n"
			"pos(x) = x > 0.0 ? x : 0.0\n"
			"theta(x) = x < 0.0 ? 0.0 : 1\n"
			"relx(x)= x< belowZeroInterval? (x-belowZeroInterval)/belowZeroInterval : (x-belowZeroInterval)/aboveZeroInterval\n\n"
			"r(x)=theta(-relx(x))*pos(-1.5*relx(x)-0.5) + theta(relx(x))*pos(2.5*relx(x)-0) + 0.9*pos(1-2.5*abs(relx(x)))\n"
			"g(x)=theta(-relx(x))*pos(-1.5*relx(x)-0.0) + theta(relx(x))*pos(2.5*relx(x)-1) + 0.9*pos(1-2.5*abs(relx(x)))\n"
			"b(x)=theta(-relx(x))*pos(-1.5*relx(x)-0.5) + theta(relx(x))*pos(2.5*relx(x)-2) + 0.9*pos(1-2.5*abs(relx(x)))\n";
		}
		else
		{
			colormap +=
					"set palette rgb 34,35,36 \n";
		}

		std::string plotscript =
		"#!/usr/bin/gnuplot\n"
		"reset\n"
		"set terminal epslatex color font 'Helvetica,10' lw 2\n"
		"set output './tmp.tex'\n\n"
		"unset key\n"
		"set grid xtics\n"
		"set xtics ("+xTics+")\n\n"
		"unset xlabel\n\n"
		"set grid ytics\n"
		"set ytics ("+yTics+")\n"
		"unset ylabel\n\n"
		+colormap+
		"set pm3d map\n"
		"sp '"+dataFile+"' binary matrix using 2:1:3 w pm3d";

		std::ofstream scriptFile( (this->get_filename()+".gp").c_str());
		if ( ! scriptFile.good() )
			error_handling::Error( std::string("Error creating file ")+(this->get_filename()+".gp"));
		scriptFile << plotscript;
		scriptFile.close();
	}
}

} /* namespace output */
} /* namespace scallop */
