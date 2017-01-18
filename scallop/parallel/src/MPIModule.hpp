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
#include <algorithm>
#include <iostream>
#include <set>

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
			{
				std::cout <<"gather: Input dimension to the target vector on the root processor must match data size!"<<std::endl;
				mpiobj.abort(2);
				std::exit(2);
			}

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
	void all_gather( T & data,  T const& dataSend )
	{
#ifdef MPI_PARALLEL
		MPIModule const& mpi = MPIModule::get_instance();
		size_t nP = mpi.get_nproc();
		auto data_ptr_send = const_cast<typename T::value_type*>(&(*dataSend.begin()));
		size_t sendcount = dataSend.end() - dataSend.begin();

		//Every processor tells how much data it intends to send
		std::vector<int> recvcounts(nP);
		recvcounts[ mpi.get_mpi_me() ] = sendcount;
		mpi.sum(recvcounts);

		//We allocate a container of the right size if data is empty
		if ( data.size() == 0 )
		{
			size_t numTotal = 0;
			for ( auto n : recvcounts )
				numTotal += n;
			data = T( numTotal );
		}

		auto data_ptr_rece = &(*data.begin());

		auto displs = recvcounts;
		displs[0] = 0;
		for ( int i = 1 ; i < static_cast<int>(nP) ; ++i )
			displs[i] = displs[i-1] + recvcounts[i-1];

		MPITypeMap<typename T::value_type> MPITypeMap;
		int ierr = MPI_Allgatherv(
				data_ptr_send,
				sendcount,
				MPITypeMap.type,
				data_ptr_rece,
				recvcounts.data(),
				displs.data(),
				MPITypeMap.type,
				MPI_COMM_WORLD);

		if ( ierr != 0 )
			std::cout << "WARNING: received exit code ierr " << ierr << std::endl;
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
		//use the single element version to check if the space is correctly allocated
		MPIModule const& mpiobj = MPIModule::get_instance();

		size_t sendcount = data.end() - data.begin();
		size_t rootSize = sendcount;
		mpiobj.bcast(rootSize,proc);
		if (sendcount != rootSize )
			data = T( rootSize );

		auto data_ptr = &(*data.begin());

		MPITypeMap<typename T::value_type> MPITypeMap;
		int ierr = MPI_Bcast(
				data_ptr,
				rootSize,
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
void MPIModule::all_to_allv( T const& send, T & recv,
		std::vector<size_t> const& sendcountsA) const
{
	MPIModule const& mpi = MPIModule::get_instance();
	size_t nP = mpi.get_nproc() ;

	std::vector<int> sendcountsi(sendcountsA.begin(),sendcountsA.end());
	auto recvcountsi = sendcountsi;
	mpi.all_to_all(recvcountsi,1);

	auto sdisplsA = sendcountsi;
	auto rdisplsA = recvcountsi;
	rdisplsA[0] = sdisplsA[0] = 0;
	for (size_t i = 1 ; i < nP; ++i )
	{
		sdisplsA[i] = sdisplsA[i-1] + sendcountsi[i-1];
		rdisplsA[i] = rdisplsA[i-1] + recvcountsi[i-1];
	}

	//We allocate a container of the right size if empty
	if ( recv.size() == 0 )
		recv = T( rdisplsA.back()+recvcountsi.back() );

	auto send_ptr = &(*send.begin());
	auto sendc_ptr = &(*sendcountsi.begin());
	auto recvc_ptr = &(*recvcountsi.begin());
	auto sdispls_ptr = &(*sdisplsA.begin());
	auto rdispls_ptr = &(*rdisplsA.begin());
	auto recv_ptr = &(*recv.begin());

	MPITypeMap< typename T::value_type > MPITypeMap;

	int ierr = MPI_Alltoallv(
			send_ptr,
			sendc_ptr,
			sdispls_ptr,
			MPITypeMap.type,
			recv_ptr,
			recvc_ptr,
			rdispls_ptr,
			MPITypeMap.type,
			MPI_COMM_WORLD);

	if ( ierr != 0 )
		std::cout << "WARNING: received exit code ierr " << ierr << std::endl;
}

template<typename T>
void MPIModule::all_to_all(T & data, size_t count) const
{
#ifdef MPI_PARALLEL
	auto data_ptr_recv = &( *data.begin() );
	MPITypeMap< typename T::value_type > MPITypeMap;

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
void MPIModule::all_gather( T & dataReceive, T const& dataSend) const
{
	delegate::all_gather_impl<T,std::is_class<T>::value > allgatherdel;
	allgatherdel.all_gather(dataReceive,dataSend);
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
void MPIModule::min(  T & data, size_t & proc ) const
{
	delegate::max_impl<T,std::is_class<T>::value > maxdel;
	T neg = T(-1.0)*data;
	maxdel.max_val( neg,proc);
	data = T(-1.0)*neg;
}

template<typename T>
void MPIModule::max(  T & data ) const
{
	size_t dummy;
	this->max( data , dummy );
}

template<typename T>
void MPIModule::min(  T & data ) const
{
	size_t dummy;
	this->min( data , dummy );
}

template<class T>
void MPIModule::save_file_parallel(
		T const& data,
		std::string const& filename,
		std::vector<size_t> ordering ) const
{
	size_t np = this->get_nproc();
	if ( ordering.empty() )
	{
		ordering = std::vector<size_t>( np );
		int i = 0 ;
		for ( auto &o : ordering )
			o = i++;
	}

	std::vector<size_t> numPerProc( np, 0 );
	numPerProc[ ordering[ this->get_mpi_me() ] ] = data.size();
	this->sum(numPerProc);
	auto displ = numPerProc;
	displ[0] = 0;
	for ( size_t ip = 1; ip < np ; ++ip )
		displ[ip] = displ[ip-1]+numPerProc[ip-1];

	MPI_Offset offset = sizeof( typename T::value_type )*displ[ ordering[ this->get_mpi_me() ] ];
	MPI_File file;
	MPI_Status status;
	MPITypeMap< typename T::value_type > MPITypeMap;

	// opening a shared file and delete if already present
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE|MPI_MODE_DELETE_ON_CLOSE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    MPI_File_close(&file);
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
	MPI_File_seek(file, offset, MPI_SEEK_SET);

	auto data_ptr = &( *data.begin() );
	MPI_File_write(file,data_ptr ,
			numPerProc[ ordering[ this->get_mpi_me() ] ],
			MPITypeMap.type,
			&status);
	MPI_File_close(&file);
}

} /* namespace parallel */
} /* namespace scallop */

#endif /* SCALLOP_PARALLEL_SRC_MPIMODULE_HPP_ */
