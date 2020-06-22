#include "Direction.hpp"
#include <stdexcept>

Direction rotate_left(Direction direction) {
	switch (direction) {
		case Direction::Right: return Direction::Up;
		case Direction::Up: return Direction::Left;
		case Direction::Left: return Direction::Down;
		case Direction::Down: return Direction::Right;
	}
	throw std::out_of_range("direction");
}
