#include "BoundingBox.h"

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
