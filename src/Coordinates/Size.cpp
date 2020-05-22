#include "Size.h"

template<class T>
T Size<T>::area() const {
	return this->width * this->height;
}

template<class T>
template<class U>
Size<T>::operator Vector2<U>() const {
	return Vector2<U>((U)this->width, (U)this->height);
}
