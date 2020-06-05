#pragma once

enum Direction {
	Right = 0,
	Up  = 1,
	Left  = 2,
	Down    = 3
};

Direction rotate_left(Direction direction);