#pragma once

#include "Vector2.h"
#include "Size.h"

template<class T>
struct BoundingBox {
	Vector2<T> start;
	Vector2<T> end;

	T width() const;
	T height() const;
	Size<T> size() const;

	bool contains(const Vector2<T> &point) const;
};

