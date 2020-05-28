#pragma once

#include "Vector2.h"
#include "Size.h"

namespace Geometry {
	template<class T>
	struct BoundingBox {
		Vector2<T> start;
		Vector2<T> end;

		T width() const;

		T height() const;

		Size<T> size() const;

		bool contains(const Vector2<T> &point) const;
	};

	template<class T>
	T BoundingBox<T>::width() const {
		return end.x - start.x;
	}

	template<class T>
	T BoundingBox<T>::height() const {
		return end.y - start.y;
	}

	template<class T>
	bool BoundingBox<T>::contains(const Vector2<T> &point) const {
		return point.x >= start.x
		       && point.x < end.x
		       && point.y >= start.y
		       && point.y < end.y;
	}

	template<class T>
	Size<T> BoundingBox<T>::size() const {
		return Size<T>(this->width(), this->height());
	}
}