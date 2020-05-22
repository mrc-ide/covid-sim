#include "Size.h"
#include <cmath>

template<class T>
T Size<T>::area() const {
	return this->width * this->height;
}

template<class T>
bool Size<T>::contains(Vector2<T> point) const {
	return point.x >= 0
		&& point.y >= 0
		&& point.x < this->width
		&& point.y < this->height;
}

template<class T>
Size<T> Size<T>::ceil() const {
	return Size<T>(std::ceil(this->width), std::ceil(this->height));
}

template<class T>
Size<T> Size<T>::operator/(Size<T> other) const {
	return Size<T>(this->width / other.width, this->height / other.height);
}

template<class T>
Size<T> Size<T>::operator*(T other) const {
	return Size<T>(this->width * other, this->height * other);
}

template<class T>
Size<T> Size<T>::operator/(T other) const {
	return Size<T>(this->width / other, this->height / other);
}

template<class T>
template<class U>
Size<T>::operator Size<U>() const {
	return Size<U>((U)this->width, (U)this->height);
}

template<class T>
template<class U>
Size<T>::operator Vector2<U>() const {
	return Vector2<U>((U)this->width, (U)this->height);
}

