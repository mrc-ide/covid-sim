#pragma once

#include "Vector2.h"

template<class T>
struct Size {
	T width;
	T height;

	Size();
	Size(T s);
	Size(T width, T height);

	T area() const;
	bool contains(const Vector2<T> &point) const;

	Size<T> floor() const;
	Size<T> ceil() const;

	Size<T> operator *(const Size<T> &other) const;
	Size<T> operator /(const Size<T> &other) const;

	Size<T> operator *(const T &other) const;
	Size<T> operator /(const T &other) const;

	template<class U>
	explicit operator Size<U>() const;
	template<class U>
	explicit operator Vector2<U>() const;
};
