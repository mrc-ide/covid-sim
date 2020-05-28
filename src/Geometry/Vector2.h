#pragma once
#include <cmath>
#include "../Dist.h"

namespace Geometry {
	template<class T>
	struct Vector2 {
		T x;
		T y;

		Vector2();

		Vector2(T x, T y);

		T length() const;

		T length_squared() const;

		T distance_to(const Vector2<T> &other) const;

		T distance_to_squared(const Vector2<T> &other) const;

		Vector2<T> abs() const;

		Vector2<T> floor() const;

		Vector2<T> ceil() const;

		Vector2<T> operator+(const Vector2<T> &other) const;

		Vector2<T> &operator+=(const Vector2<T> &other);

		Vector2<T> operator-(const Vector2<T> &other) const;

		Vector2<T> &operator-=(const Vector2<T> &other);

		Vector2<T> operator*(const Vector2<T> &other) const;

		Vector2<T> operator/(const Vector2<T> &other) const;

		Vector2<T> operator+(const T &other) const;

		Vector2<T> operator-(const T &other) const;

		Vector2<T> operator*(const T &other) const;

		Vector2<T> operator/(const T &other) const;

		template<class U>
		explicit operator Vector2<U>() const;
	};


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
	T Vector2<T>::distance_to(const Vector2<T> &other) const {
		return (T) std::sqrt(this->distance_to_squared(other));
	}

	template<class T>
	T Vector2<T>::distance_to_squared(const Vector2<T> &other) const {
		return (T) dist2_raw(this->x, this->y, other.x, other.y);
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
	Vector2<T> Vector2<T>::ceil() const {
		return Vector2<T>(std::ceil(this->x), std::ceil(this->y));
	}

	template<class T>
	Vector2<T> Vector2<T>::operator+(const Vector2<T> &other) const {
		return Vector2<T>(this->x + other.x, this->y + other.y);
	}

	template<class T>
	Vector2<T> &Vector2<T>::operator+=(const Vector2<T> &other) {
		return *this = *this + other;
	}

	template<class T>
	Vector2<T> Vector2<T>::operator-(const Vector2<T> &other) const {
		return Vector2<T>(this->x - other.x, this->y - other.y);
	}

	template<class T>
	Vector2<T> &Vector2<T>::operator-=(const Vector2<T> &other) {
		return *this = *this - other;
	}

	template<class T>
	Vector2<T> Vector2<T>::operator*(const Vector2<T> &other) const {
		return Vector2<T>(this->x * other.x, this->y * other.y);
	}

	template<class T>
	Vector2<T> Vector2<T>::operator/(const Vector2<T> &other) const {
		return Vector2<T>(this->x / other.x, this->y / other.y);
	}

	template<class T>
	Vector2<T> Vector2<T>::operator+(const T &other) const {
		return Vector2<T>(this->x + other, this->y + other);
	}

	template<class T>
	Vector2<T> Vector2<T>::operator-(const T &other) const {
		return Vector2<T>(this->x - other, this->y - other);
	}

	template<class T>
	Vector2<T> Vector2<T>::operator*(const T &other) const {
		return Vector2<T>(this->x * other, this->y * other);
	}

	template<class T>
	Vector2<T> Vector2<T>::operator/(const T &other) const {
		return Vector2<T>(this->x / other, this->y / other);
	}

	template<class T>
	template<class U>
	Vector2<T>::operator Vector2<U>() const {
		return Vector2<U>((U)
		this->x, (U)
		this->y);
	}
}