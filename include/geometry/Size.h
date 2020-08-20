#ifndef COVIDSIM_GEOMETRY_SIZE_H_INCLUDED_
#define COVIDSIM_GEOMETRY_SIZE_H_INCLUDED_

#include <cmath>
#include "Vector2.h"

namespace CovidSim { namespace Geometry {
	template<class T>
	struct Size {
		T width;
		T height;

		Size();

		explicit Size(const Vector2<T> &dimensions);

		Size(T width, T height);

		T area() const;

		bool contains(const Vector2<T> &point) const;

		Size<T> floor() const;

		Size<T> ceil() const;

		Size<T> operator*(const Size<T> &other) const;

		Size<T> operator/(const Size<T> &other) const;

		Size<T> operator*(const T &other) const;

		template<class U>
		friend Size<U> operator*(const U &other, const Size<U> &self);

		Size<T> operator/(const T &other) const;

		template<class U>
		explicit operator Size<U>() const;

		template<class U>
		explicit operator Vector2<U>() const;
	};

	template<class T>
	Size<T>::Size() : width(), height() {}

	template<class T>
	Size<T>::Size(const Vector2<T> &dimensions) : width(dimensions.x), height(dimensions.y) {}

	template<class T>
	Size<T>::Size(T _width, T _height) : width(_width), height(_height) {}

	template<class T>
	T Size<T>::area() const {
		return this->width * this->height;
	}

	template<class T>
	bool Size<T>::contains(const Vector2<T> &point) const {
		return point.x >= 0
		       && point.y >= 0
		       && point.x < this->width
		       && point.y < this->height;
	}

	template<class T>
	Size<T> Size<T>::floor() const {
		return Size<T>(std::floor(this->width), std::floor(this->height));
	}

	template<class T>
	Size<T> Size<T>::ceil() const {
		return Size<T>(std::ceil(this->width), std::ceil(this->height));
	}

	template<class T>
	Size<T> Size<T>::operator*(const Size<T> &other) const {
		return Size<T>(this->width * other.width, this->height * other.height);
	}

	template<class T>
	Size<T> Size<T>::operator/(const Size<T> &other) const {
		return Size<T>(this->width / other.width, this->height / other.height);
	}

	template<class T>
	Size<T> Size<T>::operator*(const T &other) const {
		return Size<T>(this->width * other, this->height * other);
	}

	template<class U>
	Size<U> operator*(const U &other, const Size<U> &self) {
		return Size<U>(other * self.width, other * self.height);
	}

	template<class T>
	Size<T> Size<T>::operator/(const T &other) const {
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
}}

#endif
