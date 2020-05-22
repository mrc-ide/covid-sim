#pragma once

#include "Vector2.h"

template<class T>
struct Size {
	T width;
	T height;

	T area() const;
	bool contains(Vector2<T> point) const;

	template<class U>
	explicit operator Vector2<U>() const;
};
