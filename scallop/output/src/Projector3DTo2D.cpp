/*	This file Projector3DTo2D.cpp is part of scallop.
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
 *  Created on: Mar 1, 2017
 *      Author: A. Linscheid
 */

#include "scallop/output/Projector3DTo2D.h"

namespace scallop
{
namespace output
{

void Projector3DTo2D::initialize( std::vector<size_t> grid2D,
		std::function<void(std::vector<size_t> const&,std::vector<size_t> &)> map)
{
	map_ = map;
	grid2D_ = std::move(grid2D);
	doProject_ = true;
}

std::vector<size_t> Projector3DTo2D::get_2D_grid_dimensions() const
{
	return grid2D_;
}

std::vector<size_t> const& Projector3DTo2D::operator() (std::vector<size_t> const& vect2D) const
{
	map_(vect2D,buffer3d_);
	return buffer3d_;
}

bool Projector3DTo2D::project_active() const
{
	return doProject_;
}

} /* namespace output */
} /* namespace scallop */
