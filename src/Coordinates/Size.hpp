#pragma once

template<class T>
struct Size {
	T width;
	T height;

	T area() const;
};
