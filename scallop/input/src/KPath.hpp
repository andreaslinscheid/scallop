/*	This file KPath.hpp is part of scallop.
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
 *  Created on: Nov 27, 2016
 *      Author: A. Linscheid
 */

#include "scallop/input/KPath.h"
#include <regex>

namespace scallop
{
namespace input
{

template<typename T>
void KPath<T>::read_kpath_file( std::string const& filename )
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	V kpath;
	if ( mpi.ioproc() )
	{
		std::ifstream file( filename.c_str() );
		if ( not file.good() )
			error_handling::Error(std::string()+"Problem opening kpath file "+filename,2);

		//figure out the dimension of the input k points
		std::string firstLine;
		std::getline(file,firstLine);
		std::stringstream ss(firstLine);
		std::string dummy;
		int nelem = 0;
		while( ss >> dummy)
		{
			nelem++;
		}
		dim_ = static_cast<size_t>((nelem-2-1)/2);
		if ( dim_ > 3 )
			error_handling::Error("Problem reading file: dimension figured out >3");
		file.clear();
		file.seekg(0,std::ios::beg);

		std::string label1, label2;
		std::vector<bT> v1(dim_,0);
		std::vector<bT> v2(dim_,0);
		while( true )
		{
			//read first k point
			if (not ( file >> label1 ) )
			{
				if ( file.eof() )
					break;
				scallop::error_handling::Error(std::string()+"Problem reading label1 kpath file "+filename,2);
			}

			for ( size_t i = 0 ; i < dim_; i++)
				if (not ( file >> v1[i] ) )
					scallop::error_handling::Error(std::string()+"Problem reading vector1 kpath file "+filename,2);

			size_t nkpts_sect;
			if (not ( file >> nkpts_sect ) )
				scallop::error_handling::Error(std::string()+"Problem reading section nkpts kpath file "+filename,2);

			//read second k point
			if (not ( file >> label2 ) )
				scallop::error_handling::Error(std::string()+"Problem reading label2 kpath file "+filename,2);

			for ( size_t i = 0 ; i < dim_; i++)
				if (not ( file >> v2[i] ) )
					scallop::error_handling::Error(std::string()+"Problem reading vector2 kpath file "+filename,2);

			labels_.push_back( std::make_pair(kpath.size()/dim_,label1) );
			for ( size_t iks = 0 ; iks < nkpts_sect; iks++)
				for ( size_t i = 0 ; i < dim_; i++)
				{
					bT vfbz = v1[i]+(v2[i]-v1[i])*static_cast<bT>(iks)/static_cast<bT>(nkpts_sect);
					vfbz -= std::floor(vfbz);
					kpath.push_back( vfbz );
				}
		}

		//add the last k points
		labels_.push_back( std::make_pair(kpath.size()/dim_,label2) );
		for ( size_t i = 0 ; i < dim_; i++)
		{
			bT vfbz = v2[i];
			vfbz -= std::floor(vfbz);
			kpath.push_back( vfbz );
		}
		kptsTotal_ = typename auxillary::TemplateTypedefs<T>::scallop_vector( kpath.size() );
		std::copy(kpath.begin(),kpath.end(),kptsTotal_.begin());
	}
	mpi.bcast(kptsTotal_, mpi.ioproc_index() );
	mpi.bcast(dim_, mpi.ioproc_index() );
}

template<typename T>
typename KPath<T>::V const&
KPath<T>::get_k_path() const
{
	return kptsTotal_;
}

template<typename T>
template<class DT, class OT,typename T2>
void KPath<T>::band_structure_gnuplot(
		DT const& data,
		OT const& omega,
		parallel::IrregularGridDistribution<T2> const& irrGrid,
		std::string const& filename,
		std::string quantityLabel) const
{
	assert( data.size() == omega.size()*irrGrid.get_n_pts() );

	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	//data has to be padded with the grid so that gnuplot can plot it
	std::vector< float > dataForGP( data.size()+irrGrid.get_n_pts()
							+ (mpi.get_mpi_me()==0?1+omega.size():0));

	size_t pos = 0;
	if ( mpi.get_mpi_me() == 0 )
	{
		dataForGP[pos++] = omega.size();

		for(auto w:omega)
			dataForGP[pos++] = std::real(w);
	}

	for ( size_t ik = 0; ik < irrGrid.get_n_pts(); ++ik)
	{
		dataForGP[pos++] = float(irrGrid.proc_local_to_total_index(ik))/float(irrGrid.get_n_pts_total());

		for(size_t i = 0; i < omega.size(); ++i)
			dataForGP[pos++] = -M_PI*std::imag(data[ik*omega.size()+i]);
	}

	auto dataFile = filename+".dat";
	mpi.save_file_parallel( dataForGP, dataFile);
	if ( mpi.ioproc() )
	{
		typedef typename OT::value_type aT;
		auto cmp = [](aT const& a, aT const& b){return std::real(a) < std::real(b);};
		auto mm = std::minmax_element(omega.begin(), omega.end(), cmp );
		//build the path label
		std::string xTics;
		for ( auto l : labels_)
		{
			l.second = std::regex_replace(l.second, std::regex("\\\\"), "\\\\");
			xTics += std::string("\"") + l.second
					+"\" "+std::to_string(float(l.first)/float(irrGrid.get_n_pts_total()))+",";
		}
		xTics.pop_back();

		quantityLabel = std::regex_replace(quantityLabel, std::regex("\\\\"), "\\\\");

		std::string plotscript =
		"#!/usr/bin/gnuplot\n"
		"reset\n"
		"set terminal epslatex color font 'Helvetica,10' lw 2\n"
		"set output './tmp.tex'\n\n"
		"unset key\n"
		"set grid xtics\n"
		"set xtics ("+xTics+")\n\n"
		"set ylabel \""+quantityLabel+"\" offset 0,0,0\n"
		"unset xlabel\n\n"
		"set yrange ["+std::to_string(std::real(*(mm.first)))+":"+std::to_string(std::real(*(mm.second)))+"]\n\n"
		"set arrow from 0,0 to 1,0 lw 1 lt 2 lc rgb \"black\" nohead front\n"
		"unset colorbox\n"
		"set pm3d map\n"
		"set palette model RGB functions r(gray), g(gray), b(gray)\n"
		"cutOffset(x) = x > 0.0 ? x : 0.0\n"
		"r(x)=cutOffset(-2.5*x+0)   + 0.8*cutOffset(-1+2.5*abs(x))\n"
		"g(x)=cutOffset(-2.5*x+1)   + 0.8*cutOffset(-1+2.5*abs(x))\n"
		"b(x)=cutOffset(-2.5*x+2)   + 0.8*cutOffset(-1+2.5*abs(x))\n"
		"sp '"+dataFile+"' binary matrix using 2:1:(abs($3)>exp(0)?0:-log(abs($3))>10?10:-log(abs($3))) w pm3d";

		std::ofstream scriptFile( (filename+".gp").c_str());
		if ( ! scriptFile.good() )
			error_handling::Error( std::string("Error creating file ")+(filename+".gp"));
		scriptFile << plotscript;
		scriptFile.close();
	}
}

template<typename T>
size_t
KPath<T>::get_dim() const
{
	return dim_;
}

} /* namespace input */
} /* namespace scallop */
