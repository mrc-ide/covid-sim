#ifndef COVIDSIM_GEOMETRY_VECTOR2_H_INCLUDED_
#define COVIDSIM_GEOMETRY_VECTOR2_H_INCLUDED_


#include <cmath>

namespace CovidSim { namespace Geometry {
	template<class T>
	struct DiagonalMatrix2 {
		T x;
		T y;

		DiagonalMatrix2() : x(), y() {}

		DiagonalMatrix2(T _x, T _y) : x(_x), y(_y) {}

		template<class U>
		explicit operator DiagonalMatrix2<U>() const {
			return DiagonalMatrix2<U>((U)this->x, (U)this->y);
		}
	};

	template<class T>
	struct Vector2 {
		T x;
		T y;

		Vector2();

		Vector2(T _x, T _y);

		T length() const;

		T length_squared() const;

		Vector2<T> abs() const;

		Vector2<T> operator+(const Vector2<T> &other) const;

		Vector2<T> &operator+=(const Vector2<T> &other);

		Vector2<T> operator-(const Vector2<T> &other) const;

		Vector2<T> &operator-=(const Vector2<T> &other);

		Vector2<T> operator*(const DiagonalMatrix2<T> &other) const;

		Vector2<T> operator/(const DiagonalMatrix2<T> &other) const;

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
	Vector2<T>::Vector2(T _x, T _y) : x(_x), y(_y) {}

	template<class T>
	T Vector2<T>::length() const {
		return sqrt(this->length_squared());
	}

	template<class T>
	T Vector2<T>::length_squared() const {
		return this->x * this->x + this->y * this->y;
	}

	template<class T>
	Vector2<T> Vector2<T>::abs() const {
		return Vector2<T>(std::abs(this->x), std::abs(this->y));
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
	Vector2<T> Vector2<T>::operator*(const DiagonalMatrix2<T> &other) const {
		return Vector2<T>(this->x * other.x, this->y * other.y);
	}

	template<class T>
	Vector2<T> Vector2<T>::operator/(const DiagonalMatrix2<T> &other) const {
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
		return Vector2<U>((U)this->x, (U)this->y);
	}

	typedef DiagonalMatrix2<double> DiagonalMatrix2d;
	typedef DiagonalMatrix2<float>  DiagonalMatrix2f;
	typedef DiagonalMatrix2<int>    DiagonalMatrix2i;

	typedef Vector2<double> Vector2d;
	typedef Vector2<float>  Vector2f;
	typedef Vector2<int>    Vector2i;


	Vector2d operator*(const DiagonalMatrix2d &left, const Vector2f &right);
	Vector2d operator*(const Vector2f &left, const DiagonalMatrix2d &right);

	Vector2d operator-(const Vector2d &left, const Vector2i &right);
	Vector2d operator-(const Vector2i &left, const Vector2d &right);
}}

#endif
