#pragma once

#include "Direction.h"
#include "Vector2.h"

struct MicroCellPosition : Vector2<int> {
	MicroCellPosition(int x, int y) : Vector2<int>(x, y){};

	MicroCellPosition operator+(Direction direction) const;
	MicroCellPosition& operator+=(Direction direction);
};

