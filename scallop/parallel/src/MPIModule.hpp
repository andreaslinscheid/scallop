/*	This file MPIModule.hpp is part of scallop.
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

#ifndef SCALLOP_PARALLEL_SRC_MPIMODULE_HPP_
#define SCALLOP_PARALLEL_SRC_MPIMODULE_HPP_

#include "scallop/parallel/MPIModule.h"
#include <vector>
#include <complex>
#include <type_traits>
#include <scallop/error_handling/Error.h>
#include <algorithm>

namespace scallop
{
namespace parallel
{

namespace delegate
{

template<class T, bool isclass>
struct sum_single_impl
{
};

//version for random access iterator based containers
template<class T>
struct sum_single_impl<T, true >
{
	void sum( T & data , size_t proc )
	{
#ifdef MPI_PARALLEL
		//make sure we have a random access iterator, otherwise this will not compile
		//	and IT SHOULD NOT
		size_t size = data.end() - data.begin();
		if ( size == 0)
			return;

		auto data_ptr = &(*data.begin());
		MPITypeMap<typename T::value_type> MPITypeMap;
		MPIModule const& mpiobj = MPIModule::get_instance();
		int ierr;
		if ( mpiobj.get_mpi_me() == proc)
		{
			ierr = MPI_Reduce(
					MPI_IN_PLACE,
					data_ptr,
					size,
					MPITypeMap.type,
					MPI_SUM,
					static_cast<int>(proc),
					MPI_COMM_WORLD);
		}
		else
		{
			ierr = MPI_Reduce(
					data_ptr,
					data_ptr,
					size,
					MPITypeMap.type,
					MPI_SUM,
					static_cast<int>(proc),
					MPI_COMM_WORLD);
		}

		if ( ierr != 0 )
			std::cout << "WARNING: received exit code ierr " << ierr << std::endl;
#endif
	}
};

template<class T, bool isclass>
struct sum_impl
{
};

//version for random access iterator based containers
template<class T>
struct sum_impl<T, true >
{
	void sum( T & data )
	{
#ifdef MPI_PARALLEL
		//make sure we have a random access iterator, otherwise this will not compile
		//	and IT SHOULD NOT
		size_t size = data.end() - data.begin();
		if ( size == 0)
			return;

		auto data_ptr = &(*data.begin());
		MPITypeMap<typename T::value_type> MPITypeMap;
		int ierr = MPI_Allreduce(
				MPI_IN_PLACE,
				data_ptr,
				size,
				MPITypeMap.type,
				MPI_SUM,
				MPI_COMM_WORLD);

		if ( ierr != 0 )
			std::cout << "WARNING: received exit code ierr " << ierr << std::endl;
#endif
	}
};

//version for single complex values
template<class T>
struct sum_impl< std::complex<T>, true >
{
	void sum( std::complex<T> & data )
	{
#ifdef MPI_PARALLEL
		MPITypeMap<T> MPITypeMap;
		int ierr = MPI_Allreduce(
				MPI_IN_PLACE,
				&data,
				1,
				MPITypeMap.type,
				MPI_SUM,
				MPI_COMM_WORLD);

		if ( ierr != 0 )
			std::cout << "WARNING: received exit code ierr " << ierr << std::endl;
#endif
	}
};

//version for single values
template<class T>
struct sum_impl<T, false >
{
	void sum( T & data )
	{
#ifdef MPI_PARALLEL
		MPITypeMap<T> MPITypeMap;
		int ierr = MPI_Allreduce(
				MPI_IN_PLACE,
				&data,
				1,
				MPITypeMap.type,
				MPI_SUM,
				MPI_COMM_WORLD);

		if ( ierr != 0 )
			std::cout << "WARNING: received exit code ierr " << ierr << std::endl;
#endif
	}
};

template<class T, bool isclass>
struct gather_impl
{
};

//version for random access iterator based containers
template<class T>
struct gather_impl<T, true >
{
	void gather( T & dataRoot,size_t blocksize,  T const& dataAll, size_t proc )
	{
#ifdef MPI_PARALLEL
//		//We make an assumption about dataRoot being allocated to the
//		//	correct total size. Then we assume that data is distributed
//		//	among processor equally with an overhead of N elements
//		//	placed in the first N procs
//		MPIModule const& mpiobj = MPIModule::get_instance();
//		int * recvcounts = NULL;
//		int * displs = NULL;
//
//		//Yes, the MPI library promises to NOT modify 'data_ptr_send'
//		auto data_ptr_send = const_cast<typename T::value_type*>(&(*dataAll.begin()));
//		decltype(data_ptr_send) data_ptr_recv = NULL;;
//		if ( mpiobj.get_mpi_me() == proc)
//			data_ptr_recv = &(*dataRoot.begin());
//
//		if ( mpiobj.get_mpi_me() == proc)
//		{
//			size_t sizeGather = dataRoot.end() - dataRoot.begin();
//			if ( sizeGather == 0 || sizeGather % blocksize != 0 )
//				scallop::error_handling::Error("Problem: container not allocated correct size",5);
//			std::vector< std::pair<size_t,size_t> > distr = mpiobj.distribute(sizeGather/blocksize);
//			recvcounts = new int [distr.size()];
//			displs = new int [distr.size()];
//			for ( size_t i = 0 ; i < distr.size() ; ++i )
//			{
//				recvcounts[i] = static_cast<size_t>(static_cast<int>(distr[i].second)
//						-static_cast<int>(distr[i].first))*blocksize;
//				displs[i] = distr[i].first*blocksize;
//			}
//		}
//		size_t sizeSend = dataAll.end() - dataAll.begin();
//
//		MPITypeMap<typename T::value_type> MPITypeMap;
//		int ierr = MPI_Gatherv(
//				data_ptr_send,
//				sizeSend,
//				MPITypeMap.type,
//				data_ptr_recv,
//				recvcounts,
//				displs,
//				MPITypeMap.type,
//				proc,
//				MPI_COMM_WORLD);
//
//		if ( ierr != 0 )
//			std::cout << "WARNING: received exit code ierr " << ierr << std::endl;
//
//		delete [] recvcounts;
//		delete [] displs;
#endif
	}

	void gather( T & dataRoot, T const& dataAll, size_t proc )
	{
#ifdef MPI_PARALLEL
		//We make an assumption about dataRoot being allocated to the
		//	correct total size. Then we assume that data is distributed
		//	among processor equally with an overhead of N elements
		//	placed in the first N procs
		MPIModule const& mpiobj = MPIModule::get_instance();
		int * recvcounts = NULL;
		int * displs = NULL;

		std::vector<size_t> dataProc( mpiobj.get_nproc() , 0 );
		size_t sizeSend = dataAll.end() - dataAll.begin();
		dataProc[mpiobj.get_mpi_me()] = sizeSend;
		mpiobj.sum( dataProc , proc );

		if ( mpiobj.get_mpi_me() == proc)
		{
			size_t sizeGather = dataRoot.end() - dataRoot.begin();
			size_t sizeactual = 0;
			for ( auto d : dataProc )
				sizeactual += d;
			if ( sizeGather != sizeactual )
				scallop::error_handling::Error(
						"gather: Input dimension to the target vector on the root processor must match data size!",2);

			recvcounts = new int [ mpiobj.get_nproc() ];
			displs = new int [ mpiobj.get_nproc() ];
			for ( size_t i = 0 ; i < mpiobj.get_nproc() ; ++i )
			{
				recvcounts[i] = static_cast<int>(dataProc[i]);
				displs[i] = (i==0 ? 0 : displs[i-1]+recvcounts[i-1]);
			}
		}

		//Yes, the MPI library promises to NOT modify 'data_ptr_send'
		auto data_ptr_send = const_cast<typename T::value_type*>(dataAll.data());
		decltype(data_ptr_send) data_ptr_recv = NULL;
		if ( mpiobj.get_mpi_me() == proc)
			data_ptr_recv = dataRoot.data();

		MPITypeMap<typename T::value_type> MPITypeMap;
		int ierr = MPI_Gatherv(
				data_ptr_send,
				sizeSend,
				MPITypeMap.type,
				data_ptr_recv,
				recvcounts,
				displs,
				MPITypeMap.type,
				proc,
				MPI_COMM_WORLD);

		if ( ierr != 0 )
			std::cout << "WARNING: received exit code ierr " << ierr << std::endl;

		delete [] recvcounts;
		delete [] displs;
#endif
	}
};

//version for single values
template<class T>
struct gather_impl<T, false >
{
	void gather( T & data )
	{
#ifdef MPI_PARALLEL
		MPITypeMap<typename T::value_type> MPITypeMap;
		int ierr = MPI_Allreduce(
				&data,
				&data,
				1,
				MPITypeMap.type,
				MPI_IN_PLACE,
				MPI_COMM_WORLD);

		if ( ierr != 0 )
			std::cout << "WARNING: received exit code ierr " << ierr << std::endl;
#endif
	}
};

template<class T,bool isclass>
struct scatter_impl
{
};


template<class T>
struct scatter_impl<T,true>
{
	void scatter( T const& dataRoot,size_t blocksize,  T & dataAll, size_t proc )
	{
#ifdef MPI_PARALLEL
//		//We make an assumption about dataRoot being allocated to the
//		//	correct total size. Then we assume that data is distributed
//		//	among processor equally with an overhead of N elements
//		//	placed in the first N procs
//		MPIModule const& mpiobj = MPIModule::get_instance();
//		int * sendcounts = NULL;
//		int * displs = NULL;
//
//		std::vector< size_t > sizesProc ( mpiobj.get_nproc() );
//		if ( mpiobj.get_mpi_me() == proc)
//		{
//			size_t sizeSend = dataRoot.end() - dataRoot.begin();
//			if ( sizeSend == 0 || sizeSend % blocksize != 0 )
//				scallop::error_handling::Error("Problem: container not allocated correct size",5);
//			std::vector< std::pair<size_t,size_t> > distr = mpiobj.distribute(sizeSend/blocksize);
//			sendcounts = new int [distr.size()];
//			displs = new int [distr.size()];
//			for ( size_t i = 0 ; i < distr.size() ; ++i )
//			{
//				sizesProc[i] = static_cast<size_t>(static_cast<int>(distr[i].second)
//						-static_cast<int>(distr[i].first))*blocksize;
//				sendcounts[i] = sizesProc[i];
//				displs[i] = distr[i].first*blocksize;
//			}
//		}
//
//		mpiobj.bcast(sizesProc);
//		size_t sizeRecei = dataAll.end() - dataAll.begin();
//		size_t sizeActual =  sizesProc[ mpiobj.get_mpi_me() ];
//		if ( sizeRecei != sizeActual )
//			scallop::error_handling::Error("Scatter: Receive container not allocated to correct size",100+mpiobj.get_mpi_me());
//
//		auto data_ptr_recv = &(*dataAll.begin());
//		decltype(data_ptr_recv) data_ptr_send = NULL;;
//		//Yes, the MPI library promises to NOT modify 'data_ptr_send'
//		if ( mpiobj.get_mpi_me() == proc)
//			data_ptr_send = const_cast<typename T::value_type*>(&(*dataRoot.begin()));
//
//		MPITypeMap<typename T::value_type> MPITypeMap;
//		int ierr = MPI_Scatterv(
//				data_ptr_send,
//				sendcounts,
//				displs,
//				MPITypeMap.type,
//				data_ptr_recv,
//				sizeRecei,
//				MPITypeMap.type,
//				proc,
//				MPI_COMM_WORLD);
//
//		if ( ierr != 0 )
//			std::cout << "WARNING: received exit code ierr " << ierr << std::endl;
//
//		delete [] sendcounts;
//		delete [] displs;
#endif
	}
};

template<class T,bool isclass>
struct all_gather_impl
{
};


//version for random access iterator based containers
template<class T>
struct all_gather_impl<T, true >
{
	void all_gather( T & data,  T const& dataSend, size_t eleReceive, size_t blocksize )
	{
#ifdef MPI_PARALLEL
//	This was specific to the parallization in the prelimnary version of the code. Change or remove it here.
//		//We make an assumption about data being allocated to the
//		//	correct total size on each proc. Then we assume that data is distributed
//		//	among processor equally with an overhead of N elements
//		//	placed in the first N procs
//		MPIModule const& mpiobj = MPIModule::get_instance();
//		int * recvcounts = NULL;
//		int * displs = NULL;
//
//		auto data_ptr_rece = &(*data.begin());
//		auto data_ptr_send = const_cast<typename T::value_type*>(&(*dataSend.begin()));
//		size_t sendcount = dataSend.end() - dataSend.begin();
//
//		if ( eleReceive == 0 || eleReceive % blocksize != 0 )
//			scallop::error_handling::Error("Problem: container not allocated correct size",5);
//
//		std::vector< std::pair<size_t,size_t> > distr = mpiobj.distribute(eleReceive/blocksize);
//		recvcounts = new int [distr.size()];
//		displs = new int [distr.size()];
//		for ( size_t i = 0 ; i < distr.size() ; ++i )
//		{
//			recvcounts[i] = static_cast<size_t>(static_cast<int>(distr[i].second)
//					-static_cast<int>(distr[i].first))*blocksize;
//			displs[i] = distr[i].first*blocksize;
//		}
//
//		MPITypeMap<typename T::value_type> MPITypeMap;
//		int ierr = MPI_Allgatherv(
//				data_ptr_send,
//				sendcount,
//				MPITypeMap.type,
//				data_ptr_rece,
//				recvcounts,
//				displs,
//				MPITypeMap.type,
//				MPI_COMM_WORLD);
//
//		if ( ierr != 0 )
//			std::cout << "WARNING: received exit code ierr " << ierr << std::endl;
//
//		delete [] recvcounts;
//		delete [] displs;
#endif
	}
};


template<class T,bool isclass>
struct broad_cast_impl
{
};

//version for random access iterator based containers
template<class T>
struct broad_cast_impl<T, true >
{
	void broad_cast( T & data, size_t proc )
	{
#ifdef MPI_PARALLEL
		auto data_ptr = &(*data.begin());
		size_t sendcount = data.end() - data.begin();

#ifdef DEBUG_BUILD
		//use the single element version to check if the space is correctly allocated
		MPIModule const& mpiobj = MPIModule::get_instance();
		size_t rootSize = sendcount;
		mpiobj.bcast(rootSize,proc);
		if ( sendcount != rootSize )
			scallop::error_handling::Error("broadcast called with inconsistently allocated container!",5);
#endif

		MPITypeMap<typename T::value_type> MPITypeMap;
		int ierr = MPI_Bcast(
				data_ptr,
				sendcount,
				MPITypeMap.type,
				proc,
				MPI_COMM_WORLD);

		if ( ierr != 0 )
			std::cout << "WARNING: received exit code ierr " << ierr << std::endl;
#endif
	}
};

//version for single values
template<class T>
struct broad_cast_impl<T, false >
{
	void broad_cast( T & data, size_t proc )
	{
#ifdef MPI_PARALLEL
		MPITypeMap<T> MPITypeMap;
		int ierr = MPI_Bcast(
				&data,
				1,
				MPITypeMap.type,
				proc,
				MPI_COMM_WORLD);

		if ( ierr != 0 )
			std::cout << "WARNING: received exit code ierr " << ierr << std::endl;
#endif
	}
};


template<class T,bool isclass>
struct max_impl
{
};

//version for single values
template<class T>
struct max_impl<T, false >
{
	void max_val(T & data, size_t & proc )
	{
#ifdef MPI_PARALLEL
		MPIModule const& mpi = MPIModule::get_instance();
		MPITypeMap<T> MPITypeMap;

		std::vector<T> comp( mpi.get_nproc() , T(0) );
		comp[ mpi.get_mpi_me() ] = data;
		mpi.sum( comp );
		auto itM = std::max_element( comp.begin(), comp.end() );
		proc = itM - comp.begin();
		data = *itM;
#else
		proc = 0;
#endif
	}
};

};

template<typename T>
void MPIModule::sum( T & data , size_t proc) const
{
	delegate::sum_single_impl<T,std::is_class<T>::value > sumsdel;
	sumsdel.sum( data, proc );
}

template<typename T>
void MPIModule::sum( T & data  ) const
{
	delegate::sum_impl<T,std::is_class<T>::value > sumdel;
	sumdel.sum( data );
}

template<typename T>
void MPIModule::all_to_all( typename auxillary::TemplateTypedefs<T>::scallop_vector & data, size_t count) const
{
#ifdef MPI_PARALLEL
	auto data_ptr_recv = &( *data.begin() );
	MPITypeMap<T> MPITypeMap;

	int ierr = MPI_Alltoall(
			MPI_IN_PLACE,
			static_cast<int>(count),
			MPITypeMap.type,
			data_ptr_recv,
			static_cast<int>(count),
			MPITypeMap.type,
			MPI_COMM_WORLD);

	if ( ierr != 0 )
		std::cout << "WARNING: received exit code ierr " << ierr << std::endl;
#endif
}

template<typename T>
void MPIModule::gather( T & dataRoot, size_t blocksize, T const& dataAll,  size_t proc) const
{
	delegate::gather_impl<T,std::is_class<T>::value > gatherdel;
	gatherdel.gather(dataRoot,blocksize,dataAll,proc);
}

template<typename T>
void MPIModule::gather(
		T & dataRoot, T const& dataAll,
		size_t proc ) const
{
	delegate::gather_impl<T,std::is_class<T>::value > gatherdel;
	gatherdel.gather(dataRoot,dataAll,proc);
}

template<typename T>
void MPIModule::all_gather( T & dataReceive, T const& dataSend, size_t blocksize) const
{
	this->all_gather(dataReceive,dataSend,dataReceive.size(),blocksize);
}

template<typename T>
void MPIModule::all_gather( T & dataReceive, T const& dataSend, size_t eleReceive, size_t blocksize) const
{
	delegate::all_gather_impl<T,std::is_class<T>::value > allgatherdel;
	allgatherdel.all_gather(dataReceive,dataSend,eleReceive,blocksize);
}

template<typename T>
void MPIModule::scatter( T const& dataRoot, size_t blocksize, T & dataAll,  size_t proc) const
{
	delegate::scatter_impl<T,std::is_class<T>::value > scatterdel;
	scatterdel.scatter(dataRoot,blocksize,dataAll,proc);
}

template<typename T>
void MPIModule::bcast(  T & data, size_t proc ) const
{
	delegate::broad_cast_impl<T,std::is_class<T>::value > bcastdel;
	bcastdel.broad_cast(data,proc);
}

template<typename T>
void MPIModule::max(  T & data, size_t & proc ) const
{
	delegate::max_impl<T,std::is_class<T>::value > maxdel;
	maxdel.max_val(data,proc);
}

template<typename T>
void MPIModule::max(  T & data ) const
{
	size_t dummy;
	this->max( data , dummy );
}

} /* namespace parallel */
} /* namespace scallop */

#endif /* SCALLOP_PARALLEL_SRC_MPIMODULE_HPP_ */
