/*	This file Projector3DTo2D.h is part of scallop.
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

#ifndef SCALLOP_OUTPUT_PROJECTOR3DTO2D_H_
#define SCALLOP_OUTPUT_PROJECTOR3DTO2D_H_

#include <vector>
#include <functional>

namespace scallop
{
namespace output
{

class Projector3DTo2D
{
public:
	void initialize( std::vector<size_t> grid2D,
			std::function<void(std::vector<size_t> const&,std::vector<size_t> &)> map);

	std::vector<size_t> get_2D_grid_dimensions() const;

	std::vector<size_t> const& operator() (std::vector<size_t> const& vect2D) const;

	bool project_active() const;
private:
	bool doProject_ = false;

	mutable std::vector<size_t> buffer3d_ = {0,0,0};

	std::vector<size_t> grid2D_;

	std::function<void(std::vector<size_t> const&,std::vector<size_t> &)> map_;

};

} /* namespace output */
} /* namespace scallop */

#endif /* SCALLOP_OUTPUT_PROJECTOR3DTO2D_H_ */
