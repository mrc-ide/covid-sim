#include "Vector2.h"
#include <cmath>

template<class T>
Vector2<T>::Vector2() : x(), y() {}

template<class T>
Vector2<T>::Vector2(T x, T y) : x(x), y(y) {}

template<class T>
T Vector2<T>::length() const {
	return sqrt(this->length());
}

template<class T>
T Vector2<T>::length_squared() const {
	return this->x * this->x + this->y * this->y;
}

template<class T>
T Vector2<T>::distance_to(Vector2<T> other) const {
	return (T)std::sqrt(this->distance_squared_to(other));
}

template<class T>
T Vector2<T>::distance_squared_to(Vector2<T> other) const {
	return (T)dist2_raw(this->x, this->y, other.x, other.y);
}

template<class T>
Vector2<T> Vector2<T>::abs() const {
	return Vector2<T>(std::abs(this->x), std::abs(this->y));
}

template<class T>
Vector2<T> Vector2<T>::floor() const {
	return Vector2<T>(std::floor(this->x), std::floor(this->y));
}

template<class T>
Vector2<T> Vector2<T>::operator+(Vector2<T> other) const {
	return Vector2<T>(this->x + other.x, this->y + other.y);
}

template<class T>
Vector2<T> Vector2<T>::operator+=(Vector2<T> other) {
	return *this = *this + other;
}

template<class T>
Vector2<T> Vector2<T>::operator-(Vector2<T> other) const {
	return Vector2<T>(this->x - other.x, this->y - other.y);
}

template<class T>
Vector2<T> Vector2<T>::operator-=(Vector2<T> other) {
	return *this = *this - other;
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
Vector2<T> Vector2<T>::operator+(T other) const {
	return Vector2<T>(this->x + other, this->y + other);
}

template<class T>
Vector2<T> Vector2<T>::operator-(T other) const {
	return Vector2<T>(this->x - other, this->y - other);
}

template<class T>
Vector2<T> Vector2<T>::operator*(T other) const {
	return Vector2<T>(this->x * other, this->y * other);
}

template<class T>
Vector2<T> Vector2<T>::operator/(T other) const {
	return Vector2<T>(this->x / other, this->y / other);
}

template<class T>
template<class U>
Vector2<T>::operator Vector2<U>() const {
	return Vector2<U>((U)this->x, (U)this->y);
}
