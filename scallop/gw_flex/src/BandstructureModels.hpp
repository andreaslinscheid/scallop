/*	This file BandstructureModels.hpp is part of scallop.
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
 *  Created on: Nov 26, 2016
 *      Author: A. Linscheid
 */

#include "scallop/gw_flex/BandstructureModels.h"
#include "scallop/input/InputFile.h"

namespace scallop
{
namespace gw_flex
{

template<typename T>
BandstructureModels<T>::~BandstructureModels()
{

}

template<typename T>
void TwoBandCosine<T>::compute_at_k( vbT kpts, size_t nkpts, v & unitary, vbT & energyEV ) const
{
	size_t block = 4*16;
	size_t nunitary =  nkpts*block;
	size_t nev =  nkpts*8;
	if ( unitary.size() != nunitary )
		unitary = v(nunitary);
	if ( energyEV.size() != nev )
		energyEV = vbT(nev);

	MemoryLayout meml;
	meml.initialize_layout_2pt_obj(2);

	size_t dim = kpts.size()/nkpts;
	for ( size_t ik = 0; ik < nkpts; ++ik)
	{
		bT kx =  kpts[ik*dim+0];
		bT ky =  kpts[ik*dim+1];
		for ( size_t ns = 0 ; ns < 4 ; ++ns)
		{
			energyEV[ik*8+meml.memory_layout_2pt_diagonal_nsc(0,ns)]
			         = (ns<2?1.0:-1.0) * ( 0.25*t1_*(std::cos(2.0*M_PI*kx)+std::cos(2.0*M_PI*ky)-2.0)+Eh_);

			energyEV[ik*8+meml.memory_layout_2pt_diagonal_nsc(1,ns)]
			         = (ns<2?1.0:-1.0) * ( 0.25*t2_*(std::cos(2.0*M_PI*kx)+std::cos(2.0*M_PI*ky)+2.0)+Ee_);
		}
		std::fill(&(unitary[ik*block]),&(unitary[ik*block])+block, T(0) );
		for ( size_t ibnd = 0; ibnd < 2 ; ++ibnd )
			for ( size_t ns = 0 ; ns < 4 ; ++ns)
				unitary[ik*block + meml.memory_layout_2pt_obj_nsc(ibnd,ns,ibnd,ns)] = 1.0;
	}
}

template<typename T>
size_t TwoBandCosine<T>::get_nOrb() const
{
	return 2;
}

template<typename T>
void TwoBandCosine<T>::load_model_parameters( std::istream & s )
{
	input::InputFile param;
	param.read_input_stream(s);

	std::string buf = param.get_input_config_value("t1");
	std::stringstream ss(buf);
	if ( ! (ss >> t1_) )
		t1_  = 0;
	ss.clear();

	buf = param.get_input_config_value("t2");
	ss.str(buf);
	if ( ! (ss >> t2_) )
		t2_  = 0;
	ss.clear();

	buf = param.get_input_config_value("Eh");
	ss.str(buf);
	if ( ! (ss >> Eh_) )
		Eh_  = 0;
	ss.clear();

	buf = param.get_input_config_value("Ee");
	ss.str(buf);
	if ( ! (ss >> Ee_) )
		Ee_  = 0;
}


template<typename T>
void OneBandCosine<T>::compute_at_k( vbT kpts, size_t nkpts, v & unitary, vbT & energyEV ) const
{
	size_t block = 16;
	size_t nunitary =  nkpts*block;
	size_t nev =  nkpts*4;
	if ( unitary.size() != nunitary )
		unitary = v(nunitary);
	if ( energyEV.size() != nev )
		energyEV = vbT(nev);

	MemoryLayout meml;
	meml.initialize_layout_2pt_obj(1);

	size_t dim = kpts.size()/nkpts;
	for ( size_t ik = 0; ik < nkpts; ++ik)
	{
		bT kx =  kpts[ik*dim+0];
		bT ky =  kpts[ik*dim+1];
		for ( size_t ns = 0 ; ns < 4 ; ++ns)
		{
			energyEV[ik*4+meml.memory_layout_2pt_diagonal_nsc(0,ns)]
			         = (ns<2?1.0:-1.0) * ( -2.0*t_*(std::cos(2.0*M_PI*kx)+std::cos(2.0*M_PI*ky))
			        		 -4.0*tp_*std::cos(2.0*M_PI*kx)*std::cos(2.0*M_PI*ky));
		}
		std::fill(&(unitary[ik*block]),&(unitary[ik*block])+block, T(0) );
		for ( size_t ns = 0 ; ns < 4 ; ++ns)
			unitary[ik*block + meml.memory_layout_2pt_obj_nsc(0,ns,0,ns)] = 1.0;
	}
}

template<typename T>
size_t OneBandCosine<T>::get_nOrb() const
{
	return 1;
}

template<typename T>
void OneBandCosine<T>::load_model_parameters( std::istream & s )
{
	input::InputFile param;
	param.read_input_stream(s);

	std::string buf = param.get_input_config_value("t");
	std::stringstream ss(buf);
	if ( ! (ss >> t_) )
		t_  = 0;
	ss.clear();

	buf = param.get_input_config_value("tp");
	ss.str(buf);
	if ( ! (ss >> tp_) )
		tp_  = 0;
}
} /* namespace gw_flex */
} /* namespace scallop */
