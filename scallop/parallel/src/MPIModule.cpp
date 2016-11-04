/*	This file MPIModule.cpp is part of scallop.
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
 *  Created on: Nov 2, 2016
 *      Author: A. Linscheid
 */

#include "scallop/parallel/MPIModule.h"
#ifdef MPI_PARALLEL
#endif

namespace scallop
{
namespace parallel
{

MPIModule& MPIModule::get_instance()
{
	static MPIModule instance;
	return instance;
}

MPIModule::MPIModule()
{
	this->mpi_init();
}

void MPIModule::barrier() const
{
#ifdef MPI_PARALLEL
	MPI_Barrier( MPI_COMM_WORLD );
#endif
}

size_t MPIModule::ioproc_index() const
{
	return 0 ;
}

bool MPIModule::ioproc () const
{
	return mpiMe_ == this->ioproc_index() ;
}

bool MPIModule::true_on_all(bool val) const
{
	int v = val ? 0 : 1;
#ifdef MPI_PARALLEL
	this->sum( v );
#endif
	return v == 0;
}

bool MPIModule::false_on_all(bool val) const
{
	int v = val ? 1 : 0;
#ifdef MPI_PARALLEL
	this->sum( v );
#endif
	return v != 0;
}

void MPIModule::mpi_init()
{
#ifdef MPI_PARALLEL
	// Initialize the MPI environment
    MPI_Init(NULL, NULL);

    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    nproc_ = static_cast<size_t>(worldSize);

    // Get the rank of the process
    int worldRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    mpiMe_ = static_cast<size_t>(worldRank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    name_ = std::string( processor_name );

    if ( mpiMe_ == 0 )
    	std::cout << "Running on " << nproc_ << " processors!" << std::endl;
#endif
}

size_t MPIModule::get_mpi_me() const
{
	return mpiMe_;
}

size_t MPIModule::get_nproc() const
{
	return nproc_;
}

std::string const& MPIModule::get_proc_name() const
{
	return name_;
}

MPIModule::~MPIModule()
{
#ifdef MPI_PARALLEL
	// Finalize the MPI environment.
	MPI_Finalize();
#endif
}

} /* namespace parallel */
} /* namespace scallop */
