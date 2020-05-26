#include "Param.h"

using namespace Geometry;

Size<int> Param::number_of_micro_cells() const {
	return this->number_of_cells * this->NMCL;
}

MicroCellPosition Param::get_micro_cell_position_from_cell_index(int cell_index) const {
	int x = cell_index / this->number_of_micro_cells().height;
	int y = cell_index % this->number_of_micro_cells().height;
	return {x, y};
}

bool Param::is_in_bounds(MicroCellPosition position) const {
	return position.x >= 0
		&& position.y >= 0
		&& position.x < this->number_of_micro_cells().width
		&& position.y < this->number_of_micro_cells().height;
}

int Param::get_micro_cell_index_from_position(MicroCellPosition position) const {
	int x = (position.x + this->number_of_micro_cells().width) % this->number_of_micro_cells().width;
	int y = (position.y + this->number_of_micro_cells().height) % this->number_of_micro_cells().height;
	return x * this->number_of_micro_cells().height + y;
}
