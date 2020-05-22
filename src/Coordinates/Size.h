#pragma once

#include "Vector2.h"

template<class T>
struct Size {
	T width;
	T height;

	T area() const;
	bool contains(Vector2<T> point) const;

	Size<T> ceil() const;

	Size<T> operator /(Size<T> other) const;

	Size<T> operator *(T other) const;
	Size<T> operator /(T other) const;

	template<class U>
	explicit operator Size<U>() const;
	template<class U>
	explicit operator Vector2<U>() const;
};
