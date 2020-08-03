#ifndef COVIDSIM_DIRECTION_HPP_INCLUDED_
#define COVIDSIM_DIRECTION_HPP_INCLUDED_

enum Direction {
	Right = 0,
	Up  = 1,
	Left  = 2,
	Down    = 3
};

Direction rotate_left(Direction direction);

#endif
