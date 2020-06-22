#pragma once

enum struct Direction {
	Right,
	Up,
	Left,
	Down
};

Direction rotate_left(Direction direction);
