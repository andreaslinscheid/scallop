/*	This file BandstructureModels.h is part of scallop.
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

#ifndef SCALLOP_GW_FLEX_BANDSTRUCTUREMODELS_H_
#define SCALLOP_GW_FLEX_BANDSTRUCTUREMODELS_H_

namespace scallop
{
namespace gw_flex
{

template<typename T>
class BandstructureModels
{
public:
	typedef typename auxillary::TypeMapComplex<T>::type bT;
	typedef typename auxillary::TemplateTypedefs<T>::scallop_vector v;
	typedef typename auxillary::TemplateTypedefs<bT>::scallop_vector vbT;

	virtual ~BandstructureModels();

	virtual void compute_at_k( vbT kpts, size_t nkpts, v & unitary, vbT & energyEV ) const = 0;

	virtual size_t get_nOrb( ) const = 0;

	virtual void load_model_parameters( std::istream & s ) = 0;
};

template<typename T>
class TwoBandCosine : public BandstructureModels<T>
{
public:
	using typename BandstructureModels<T>::bT;
	using typename BandstructureModels<T>::v;
	using typename BandstructureModels<T>::vbT;

	void compute_at_k( vbT kpts, size_t nkpts, v & unitary, vbT & energyEV ) const;

	size_t get_nOrb() const;

	void load_model_parameters( std::istream & s );

private:

	bT t1_ = bT(0);
	bT t2_ = bT(0);
	bT Eh_ = bT(0);
	bT Ee_ = bT(0);

};

template<typename T>
class OneBandCosine : public BandstructureModels<T>
{
public:
	using typename BandstructureModels<T>::bT;
	using typename BandstructureModels<T>::v;
	using typename BandstructureModels<T>::vbT;

	void compute_at_k( vbT kpts, size_t nkpts, v & unitary, vbT & energyEV ) const;

	size_t get_nOrb() const;

	void load_model_parameters( std::istream & s );

private:

	bT t_ = bT(0);
	bT tp_ = bT(0);
};

} /* namespace gw_flex */
} /* namespace scallop */

#include "scallop/gw_flex/src/BandstructureModels.hpp"
#endif /* SCALLOP_GW_FLEX_BANDSTRUCTUREMODELS_H_ */
