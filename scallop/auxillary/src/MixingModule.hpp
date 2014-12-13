/*	This file MixingModule.hpp is part of scallop.
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
 *  Created on: Nov 30, 2014
 *      Author: Andreas Linscheid
 */
#include "scallop/auxillary/MixingModule.h"
#include "scallop/error_handling/Error.h"

namespace scallop {
namespace auxillary {

template<class T>
MixingModule<T>::MixingModule(T const& quantityToMix,
		size_t numberIterationsConsidered)
	:	_history(numberIterationsConsidered,quantityToMix) {
	_locationConologicalHistory.clear();
	_locationConologicalHistory.reserve(numberIterationsConsidered);
	for ( size_t i = 0 ; i < _history.size() ; ++i)
		_locationConologicalHistory.push_back(i);
};

template<class T>
void MixingModule<T>::init(T const& quantityToMix,
		size_t numberIterationsConsidered) {
	if ( numberIterationsConsidered > 0)
		error_handling::Error("Mixing must use 1 or more previous steps",2);
	_history = std::vector<T>(numberIterationsConsidered,quantityToMix);
	_locationConologicalHistory.clear();
	_locationConologicalHistory.reserve(numberIterationsConsidered);
	for ( size_t i = 0 ; i < _history.size() ; ++i)
		_locationConologicalHistory.push_back(i);

	//save init parameters
	_maxNumberOfIterationsConsideredForMixing = numberIterationsConsidered;
	_vectorDim = vectorDimension;

	//set the weights to the values suggested in the paper
	_weightsForErrorIteration =
			std::vector<typename T::value_type>(numberIterationsConsidered,1.0);
	_weightsForErrorInverseJacobian = 0.01;

	//init storages
	_deltaFunctional = _history;
	_deltaVector = _history;

	if ( _beta != NULL)
		delete [] _beta;
	_beta = new T [_maxNumberOfIterationsConsideredForMixing*_maxNumberOfIterationsConsideredForMixing];
	if ( _cOfEq15b != NULL)
		delete [] _cOfEq15b;
	_cOfEq15b = new T [_maxNumberOfIterationsConsideredForMixing];
	//
	//note that we are initialized
	_isInitialized = true;

}

template<class T>
T const& MixingModule<T>::get_current_iteration() const {
	return this->get_privious_iteration(0);
}

template<class T>
T const& MixingModule<T>::get_privious_iteration(size_t numItInThePast=1) const {
	return _history[_locationConologicalHistory[numItInThePast]];
}



template<typename T>
void MixingModule<T>::mixing(double mixingParameter,T & newAndMixed){

	//save the vector from the previous iteration
	size_t const historyIndexToBeReplaced = _locationConologicalHistory.back();
	_history[historyIndexToBeReplaced].swap(newAndMixed);

	//cycle the indices
	for ( size_t i = _locationConologicalHistory.size()-1 ; i > 0 ; i--)
		_locationConologicalHistory[i] = _locationConologicalHistory[i-1];
	_locationConologicalHistory.front() =  historyIndexToBeReplaced;

	//determine the number of steps in history available to improve the mixing
	size_t const numberOfIterationsConsideredForMixing =
			_numberOfIterationsSoFar  < _maxNumberOfIterationsConsideredForMixing ?
					_numberOfIterationsSoFar : _maxNumberOfIterationsConsideredForMixing;

	for ( size_t i = 0 ; i < _deltaFunctional[historyIndexToBeReplaced].size(); ++i)
		_deltaFunctional[historyIndexToBeReplaced][i] =
				this->get_current_iteration()[i]-this->get_privious_iteration(1)[i];

	//linear mixing so far
	for ( size_t i = 0 ; i < _deltaFunctional[historyIndexToBeReplaced].size(); ++i)
		newAndMixed[i] = this->get_privious_iteration(1)[i] +
			mixingParameter*_deltaFunctional[historyIndexToBeReplaced][i];

	//determine the corrections to the linear mixing
	for ( size_t iPast = 0 ; iPast < numberOfIterationsConsideredForMixing; ++iPast){
		for ( size_t jPast = iPast ; jPast < numberOfIterationsConsideredForMixing; ++jPast){

			size_t aII = _locationConologicalHistory[iPast];
			size_t aIJ = _locationConologicalHistory[jPast];
			typename T::value_type dotProduct = 0;
			for ( size_t i = 0 ; i < _deltaFunctional[aII].size(); ++i)
				dotProduct += _deltaFunctional[aII][i]*_deltaFunctional[aIJ][i];

			//formula for a Eq. (13a)
			_beta[iPast*numberOfIterationsConsideredForMixing+jPast] =
					_weightsForErrorIteration[iPast]*_weightsForErrorIteration[jPast]*dotProduct;
		}
	}
	for ( size_t iPast = 0 ; iPast < numberOfIterationsConsideredForMixing; ++iPast){
		for ( size_t jPast = iPast+1 ; jPast < numberOfIterationsConsideredForMixing; ++jPast){
			_beta[jPast*numberOfIterationsConsideredForMixing+iPast] =
					_beta[iPast*numberOfIterationsConsideredForMixing+jPast];
		}
	}
//	_dotProductBuffer[historyIndexToBeReplaced*_maxNumberOfIterationsConsideredForMixing + i]
//	                  = 0;
//	for (size_t i = 0; i < _maxNumberOfIterationsConsideredForMixing; ++i)
//		_dotProductBuffer[historyIndexToBeReplaced*_maxNumberOfIterationsConsideredForMixing + i]
//		                  +=
//
//	//generate F (from the introduction text F=v1-v2)
//	// pointerToVectorOld is now F
//	for (size_t i = _vectorDim;i--;)
//		pointerToVectorOld[i] -= pointerToVectorNewAndOut[i];
//	//
//	//generate dF and dn from Eq. (9) (only used if we already have passed one iteration)
//	//See also note at the bottom of this routine to understand why the latest entry of
//	//	|dF> contains F here and |dn> contains n from the previous iteration.
//	if ( numberOfIterationsConsideredForMixing != 0){
//		for (size_t i=_vectorDim;i--;) {
//			_deltaFunctional(i,0) = pointerToVectorOld[i]/*F*/-_deltaFunctional(i,0)/*F previous*/;
//			_deltaVector(i,0) = pointerToVectorNewAndOut[i]/*n*/-_deltaVector(i,0)/* n previous*/;
//		}
//		//
//		//normalize
//		double norm = 0.0;
//		for (size_t i=_vectorDim;i--;)
//			norm += std::pow(std::abs(_deltaFunctional(i,0)),2);
//	//s		norm += std::pow(std::abs(_deltaVector(i,0)),2);
//		norm = 1.0/std::sqrt(norm);
//		//
//		for (size_t i=_vectorDim;i--;)
//			_deltaFunctional(i,0) *= norm;
//		for (size_t i=_vectorDim;i--;)
//			_deltaVector(i,0) *= norm;
//	}
//	//
//	//compute beta of Eq. (13a):
//	//	First construct the inverse of beta
//	std::MatrixTools matTools;
//	for ( size_t i=0;i<numberOfIterationsConsideredForMixing;i++){
//		//
//		//use that <dF(i)|dF(i)>=1
//		_beta[i*numberOfIterationsConsideredForMixing+i] =
//				std::pow(_weightsForErrorInverseJacobian,2)
//				+ std::pow(_weightsForErrorIteration[i],2);
//		//
//		//use the fact that beta^-1 is symmetric
//		for ( size_t j=i+1;j<numberOfIterationsConsideredForMixing;j++){
//			_beta[i*numberOfIterationsConsideredForMixing+j] =
//					//
//					//formula for a Eq. (13a)
//					_weightsForErrorIteration[i]*_weightsForErrorIteration[j]*
//					matTools.dot_product(_vectorDim,_deltaFunctional[i],_deltaFunctional[j]);
//		}
////		for ( size_t j=0;j<i;j++){
////			_beta[i*numberOfIterationsConsideredForMixing+j]=
////					_beta[j*numberOfIterationsConsideredForMixing+i];
////		}
//	}
//	if ( numberOfIterationsConsideredForMixing > 0){
//		matTools.invert_symmetric(_beta,numberOfIterationsConsideredForMixing);
//	}
//	for ( size_t i=0;i<numberOfIterationsConsideredForMixing;i++){
//		for ( size_t j=i+1;j<numberOfIterationsConsideredForMixing;j++){
//			_beta[j*numberOfIterationsConsideredForMixing+i]=
//					_beta[i*numberOfIterationsConsideredForMixing+j];
//		}
//	}
//	for ( size_t i=numberOfIterationsConsideredForMixing;i--;){
//		_cOfEq15b[i] = matTools.dot_product(_vectorDim,_deltaFunctional[i],pointerToVectorOld/*<--F*/);
//	}
//	//
//	//fist part of Eq. 15a
//	for (size_t i = _vectorDim;i--;)
//		pointerToVectorNewAndOut[i] += mixingParameter*pointerToVectorOld[i]/*alpha * F*/;
//	//
//	//second part of Eq. 15a
//	for ( size_t i=0;i<numberOfIterationsConsideredForMixing;i++){
//		//
//		//compute gamma on the fly
//		T gammaM = 0.0;
//		for ( size_t j=0;j<numberOfIterationsConsideredForMixing;j++){
//			gammaM += _beta[i*numberOfIterationsConsideredForMixing+j]*_weightsForErrorIteration[j]*_cOfEq15b[j];
//		}
//		//
//		//compute |u> on the fly according to Eq. 13b
//		for (size_t vectorIndex = _vectorDim;vectorIndex--;)
//			pointerToVectorNewAndOut[vectorIndex] -= gammaM * _weightsForErrorIteration[i] *
//							/*|u>*/(mixingParameter*_deltaFunctional(vectorIndex,i)+_deltaVector(vectorIndex,i));
//	}
//	_numberOfIterationsSoFar++;
//	//
//	//do not save data if we consider no past iterations
//	if ( _maxNumberOfIterationsConsideredForMixing == 0 )
//		return;
//	//
//	//save for the next iteration the previous vector |n>
//	_deltaVector.push(_saveOfVectorFromPreviousIteration);
//	//
//	//and the difference of current and previous vector |F>
//	_deltaFunctional.push(pointerToVectorOld);
//	//NOTE:
//	//Here _deltaVector and _deltaFunctional do not contain the what they are called like
//	//	namely not |dF> and |dn> in the notation of the paper.
//	//	Both will be normalized at the beginning of the next iteration and then restore
	//	their naming convention. While it would be cleaner to have independent places
	//	for F and n we choose to optimize this additional allocation.
}
} /* namespace auxillary */
} /* namespace scallop */
