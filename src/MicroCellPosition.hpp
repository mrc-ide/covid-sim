#ifndef COVIDSIM_MICRO_CELL_POSITION_HPP_INCLUDED_
#define COVIDSIM_MICRO_CELL_POSITION_HPP_INCLUDED_

#include "Direction.hpp"
#include "Error.h"

struct MicroCellPosition
{
	int x;
	int y;

	inline MicroCellPosition& operator+=(Direction direction)
	{
		switch (direction) {
			case Direction::Right: this->x += 1; break;
			case Direction::Up: this->y -= 1; break;
			case Direction::Left: this->x -= 1; break;
			case Direction::Down: this->y += 1; break;
			default: ERR_CRITICAL("Unknown direction");
		}
		return *this;
	}
};

#endif
