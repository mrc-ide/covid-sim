#pragma once

#include "Vector2.h"

template<class T>
struct Size {
	T width;
	T height;

	T area() const;

	template<class U>
	explicit operator Vector2<U>() const;
};
