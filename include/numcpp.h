#ifndef NUMCPP_H_
#define NUMCPP_H_

#include <complex>
#include <iostream>

template<typename T>
namespace numcpp {

    class ndarray {

        public:

	ndarray (int n1);
	ndarray (int n1, int n2);
	ndarray (int n1, int n2, int n3);
	ndarray (int n1, int n2, int n3, int n4);

        private:
	ndarray shape();

    };

    ndarray<T>& empty (int n1); 
    ndarray<T>& empty (int n1, int n2); 
    ndarray<T>& empty (int n1, int n2, int n3); 
    ndarray<T>& zeros (int n1); 
    ndarray<T>& zeros (int n1, int n2); 
    ndarray<T>& zeros (int n1, int n2, int n3); 
}

#endif
