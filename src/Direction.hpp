#ifndef COVIDSIM_DIRECTION_HPP_INCLUDED_
#define COVIDSIM_DIRECTION_HPP_INCLUDED_

enum struct Direction {
	Right,
	Up,
	Left,
	Down
};

Direction rotate_left(Direction direction);

#endif
