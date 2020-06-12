#pragma once

#include <stdexcept>

#include "Direction.hpp"

struct MicroCellPosition
{
	int x;
	int y;

	inline MicroCellPosition& operator+=(Direction direction)
	{
		switch (direction) {
			case Right: this->x += 1; break;
			case Up: this->y -= 1; break;
			case Left: this->x -= 1; break;
			case Down: this->y += 1; break;
			default: throw std::out_of_range("direction");
		}
		return *this;
	}
};
