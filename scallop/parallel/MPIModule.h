/*	This file MPIModule.h is part of scallop.
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

#ifndef SCALLOP_PARALLEL_MPIMODULE_H_
#define SCALLOP_PARALLEL_MPIMODULE_H_

#include "scallop/auxillary/TemplateTypedefs.h"
#include <mpi.h>
#include <complex>
#include <vector>

namespace scallop
{
namespace parallel
{
/**
 * 	Interface to the MPI library.
 *
 * 	Implemented as a Singleton pattern, since there must not be two
 * 	instances of this object. Ever.
 */
class MPIModule
{
public:
	static MPIModule& get_instance();

	//Singleton: don't implement
	void operator= (MPIModule const& other) = delete;
	MPIModule(MPIModule const& other) 		= delete;

	~MPIModule();

	//returns true if this processor is the one producing i/o to the terminal
	bool ioproc () const;

	//returns the processor index that produces i/o
	size_t ioproc_index () const;

	//return the index of this processors
	size_t get_mpi_me() const;

	//return total number of processors
	size_t get_nproc() const;

	std::string const& get_proc_name() const;

	/**
	 * Performs an 'in place' all-to-all communication.
	 *
	 * @param data		Data to be redistributed
	 * @param count		Each process is sending/receiving count elements
	 */
	template<typename T>
	void all_to_all( T & data, size_t count) const;

	/**
	 * Performs an all-to-all vector communication.
	 *
	 * This high-level interface does not allow to specify displacements.
	 *
	 * @param send			Data to be send
	 * @param recv			Received data goes in here. Will be allocated to the right size, if recv.size()==0.
	 * @param sendcounts	Array of size #processors. (a=)sendcounts[n] means a elements will be send to proc n.
	 */
	template<typename T>
	void all_to_allv( T const& send, T & recv,
			std::vector<size_t> const& sendcounts) const;

	//wait for everybody to get here
	void barrier() const;

	//returns true, if the variable val is true on all processors
	bool true_on_all(bool val) const;

	//returns false, if the variable val is false on all processors
	bool false_on_all(bool val) const;

	//sum up data among all processors and distribute the result into data
	template<typename T>
	void sum( T & data ) const;

	//sum up data among all processors and place the result into data at processor 'proc'
	template<typename T>
	void sum( T & data , size_t proc) const;

	template<typename T>
	void gather(  T & dataRoot, size_t blocksize, T const& dataAll, size_t proc = 0 ) const;

	//Interface that allows specify the number of elements dataAll on each processors
	// layout is only significant at the root proc
	template<typename T>
	void gather(  T & dataRoot, T const& dataAll, size_t proc = 0 ) const;

	template<typename T>
	void all_gather( T & dataReceive, T const& dataSend) const;

	//Scatter the data in dataRoot from the proc to the other processors
	template<typename T>
	void scatter(  T const& dataRoot, size_t blocksize, T & dataAll, size_t proc = 0 ) const;

	//Broadcast the data in proc to the other processors
	template<typename T>
	void bcast(  T & data, size_t proc = 0 ) const;

	//Determine the maximal value among processors and overwrite every individual value with it.
	//Sets proc to the processor that owned the maximal value
	template<typename T>
	void max(  T & data, size_t & proc ) const;

	//Determine the maximal value among processors and overwrite every individual value with it.
	template<typename T>
	void max(  T & data ) const;

	void abort(int ierr) const;

	/**
	 * Save distributed data to a file.
	 *
	 * @param data		The distributed data on each processor.
	 * @param filename	The name of the new file.
	 * @param ordering	Optional. An array of size of # processors,
	 * 					selecting the order in which processor writes its data to the file.
	 * 					By default, every processor appends the data according to its index.
	 */
	template<class T>
	void save_file_parallel(
			T const& data,
			std::string const& filename,
			std::vector<size_t> ordering = std::vector<size_t>() ) const;
private:

	//the rank in MPI language
	size_t mpiMe_ = 0;

	size_t nproc_ = 1;

	std::string name_ = "default";

	MPIModule();

	void mpi_init();
};

template<typename T>
struct MPITypeMap
{
//	default does not define a type! Specify all types explicitly!
};

#ifdef MPI_PARALLEL
template<>
struct MPITypeMap< bool >
{
	MPI_Datatype  type = MPI_BYTE;
};

template<>
struct MPITypeMap< int >
{
	MPI_Datatype  type = MPI_INT;
};

template<>
struct MPITypeMap< size_t >
{
	MPI_Datatype  type = MPI_UINT64_T;
};

template<>
struct MPITypeMap< float >
{
	MPI_Datatype  type = MPI_FLOAT;
};

template<>
struct MPITypeMap< std::complex<float> >
{
	MPI_Datatype  type = MPI_COMPLEX;
};

template<>
struct MPITypeMap< double >
{
	MPI_Datatype  type = MPI_DOUBLE;
};

template<>
struct MPITypeMap< std::complex<double> >
{
	MPI_Datatype  type = MPI_DOUBLE_COMPLEX;
};
#endif

} /* namespace parallel */
} /* namespace scallop */

#include "scallop/parallel/src/MPIModule.hpp"
#endif /* SCALLOP_PARALLEL_MPIMODULE_H_ */
