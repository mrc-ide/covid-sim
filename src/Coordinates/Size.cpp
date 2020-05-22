#include "Size.hpp"

template<class T>
T Size<T>::area() const {
	return this->width * this->height;
}