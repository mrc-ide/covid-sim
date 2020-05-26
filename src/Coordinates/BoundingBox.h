#pragma once

#include "Vector2.h"

template<class T>
struct BoundingBox {
	Vector2<T> start;
	Vector2<T> end;

	T width() const;
	T height() const;
};

