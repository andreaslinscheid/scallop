/*	This file GapFileReader.cpp is part of scallop.
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
 *  Created on: Jan 18, 2017
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/GapFileReader.h"
#include "scallop/parallel/MPIModule.h"
#include "scallop/error_handling/Error.h"
#include "scallop/output/TerminalOut.h"
#include <fstream>
#include <assert.h>

namespace scallop
{
namespace gw_flex
{

void GapFileReader::read_file(std::string const& filename )
{
	parallel::MPIModule const& mpi = parallel::MPIModule::get_instance();
	output::TerminalOut msg;

	channels_.clear();

	//format is described in the manual
	std::ifstream file;
	std::string s;
	if ( mpi.ioproc() )
	{
		file.open( filename.c_str() );
		if ( ! file.good() )
			error_handling::Error(std::string()+"Unable to open file "+filename,1);
		msg << "Reading file " << filename;

		file >> s;
		//We'll deal with the python thingy later..
		if ( s.compare("script") == 0 )
			error_handling::Error("Not yet implemented - no python support for input gap");
		file.close();
	}
	mpi.bcast(s, mpi.ioproc_index() );

	if ( mpi.ioproc() )
	{
		file.open( filename.c_str() );
		if ( ! file.good() )
			error_handling::Error(std::string()+"Unable to open file "+filename,1);
	}

	while ( true )
	{
		bool read_file = true;
		if ( mpi.ioproc() )
			if ( ! std::getline(file,s) )
				read_file = false;

		mpi.bcast(read_file, mpi.ioproc_index() );
		mpi.bcast(s, mpi.ioproc_index() );

		if ( not read_file )
			break;

		ChannelFuncPara c;
		std::stringstream ss(s);

		std::string spinState;
		ss >> spinState;
		if ( spinState.compare("s") == 0 )
			c.spinState_ = 0;
		if ( spinState.compare("tx") == 0 )
			c.spinState_ = 1;
		if ( spinState.compare("ty") == 0 )
			c.spinState_ = 2;
		if ( spinState.compare("tz") == 0 )
			c.spinState_ = 3;

		std::string channel;
		ss >> channel;
		if ( channel.compare("s") == 0 )
			c.shape_ = &GapFileReader::s_wave;

		if ( channel.compare("px") == 0 )
			c.shape_ = &GapFileReader::px_wave;
		if ( channel.compare("py") == 0 )
			c.shape_ = &GapFileReader::py_wave;
		if ( channel.compare("pz") == 0 )
			c.shape_ = &GapFileReader::pz_wave;

		if ( channel.compare("dxy") == 0 )
			c.shape_ = &GapFileReader::dxy_wave;
		if ( channel.compare("dyz") == 0 )
			c.shape_ = &GapFileReader::dyz_wave;
		if ( channel.compare("dz2") == 0 )
			c.shape_ = &GapFileReader::dz2_wave;
		if ( channel.compare("dxz") == 0 )
			c.shape_ = &GapFileReader::dxz_wave;
		if ( channel.compare("dx2my2") == 0 )
			c.shape_ = &GapFileReader::dx2my2_wave;

		ss >> c.width_;

		double p;
		while ( ss >> p )
		{
			c.orbFactors_.push_back(p);
		}

		channels_.push_back( std::move(c) );
	}
};

double GapFileReader::radius( std::vector<double> const& k ) const
{

	return std::sqrt(this->radius2(k));
}

double GapFileReader::radius2( std::vector<double> const& k ) const
{
	double r2 = 0;
	for ( double kxi : k )
		r2 += std::pow(kxi,2);
	if ( r2 == 0 )
		r2 = 1e-12;
	return r2;
}

double GapFileReader::s_wave( std::vector<double> const& k ) const
{
	return 0.5*std::sqrt(1.0/M_PI);
}

double GapFileReader::px_wave( std::vector<double> const& k ) const
{

	return std::sqrt(3.0/(4.0*M_PI))*k[0]/this->radius(k);
}

double GapFileReader::py_wave( std::vector<double> const& k ) const
{
	assert( k.size() > 1 );
	return std::sqrt(3.0/(4.0*M_PI))*k[1]/this->radius(k);
}

double GapFileReader::pz_wave( std::vector<double> const& k ) const
{
	assert( k.size() > 2 );
	return std::sqrt(3.0/(4.0*M_PI))*k[2]/this->radius(k);
}

double GapFileReader::dxy_wave( std::vector<double> const& k ) const
{
	assert( k.size() > 1 );
	return std::sqrt(15.0/(4.0*M_PI))*k[0]*k[1]/this->radius2(k);
}

double GapFileReader::dyz_wave( std::vector<double> const& k ) const
{
	assert( k.size() > 2 );
	return std::sqrt(15.0/(4.0*M_PI))*k[2]*k[1]/this->radius2(k);
}

double GapFileReader::dz2_wave( std::vector<double> const& k ) const
{
	assert( k.size() > 2 );
	return std::sqrt(5.0/(16.0*M_PI))*(2.0*k[2]*k[2]-k[1]*k[1]-k[0]*k[0])/this->radius2(k);
}

double GapFileReader::dxz_wave( std::vector<double> const& k ) const
{
	assert( k.size() > 2 );
	return std::sqrt(15.0/(4.0*M_PI))*k[2]*k[0]/this->radius2(k);
}

double GapFileReader::dx2my2_wave( std::vector<double> const& k ) const
{
	assert( k.size() > 1 );
	return std::sqrt(15.0/(16.0*M_PI))*(k[0]*k[0]-k[1]*k[1])/this->radius2(k);
}

} /* namespace gw_flex */
} /* namespace scallop */
