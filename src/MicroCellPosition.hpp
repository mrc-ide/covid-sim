#pragma once

#include "Direction.hpp"

class MicroCellPosition {
public:
	int x;
	int y;

	MicroCellPosition operator+(Direction direction) const;
	MicroCellPosition& operator+=(Direction direction);
};

