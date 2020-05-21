#pragma once

#include "Direction.hpp"
#include "Vector2.hpp"

struct MicroCellPosition : Vector2<int> {
	MicroCellPosition(int x, int y) : Vector2<int>(x, y){};

	MicroCellPosition operator+(Direction direction) const;
	MicroCellPosition& operator+=(Direction direction);
};

