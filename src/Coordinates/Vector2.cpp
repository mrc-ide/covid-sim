#include "Vector2.hpp"
#include <cmath>

template<class T>
Vector2<T>::Vector2(T x, T y) :
	x(x),
	y(y)
{}

template<class T>
Vector2<T> Vector2<T>::abs() const {
	return Vector2<T>(std::abs(this->x), std::abs(this->y));
}

template<class T>
Vector2<T> Vector2<T>::operator+(Vector2<T> other) const {
	return Vector2<T>(this->x + other.x, this->y + other.y);
}

template<class T>
Vector2<T> Vector2<T>::operator-(Vector2<T> other) const {
	return Vector2<T>(this->x - other.x, this->y - other.y);
}

template<class T>
Vector2<T> Vector2<T>::operator*(Vector2<T> other) const {
	return Vector2<T>(this->x * other.x, this->y * other.y);
}

template<class T>
Vector2<T> Vector2<T>::operator/(Vector2<T> other) const {
	return Vector2<T>(this->x / other.x, this->y / other.y);
}

template<class T>
template<class U>
Vector2<T>::operator Vector2<U>() const {
	return Vector2<U>((U)this->x, (U)this->y);
}
