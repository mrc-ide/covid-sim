#include "BoundingBox.h"

template<class T>
T BoundingBox<T>::width() const {
	return end.x - start.x;
}

template<class T>
T BoundingBox<T>::height() const {
	return end.y - start.y;
}