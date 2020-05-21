#include "Direction.hpp"
#include <stdexcept>

Direction rotate_left(Direction direction) {
	switch (direction) {
		case Right: return Up;
		case Up: return Left;
		case Left: return Down;
		case Down: return Right;
	}
	throw std::out_of_range("direction");
}
