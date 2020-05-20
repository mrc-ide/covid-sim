#include "Param.h"

int Param::get_number_of_micro_cells_wide() const {
	return this->ncw * this->NMCL;
}

int Param::get_number_of_micro_cells_high() const {
	return this->nch * this->NMCL;
}

MicroCellPosition Param::get_micro_cell_position_from_cell_index(int cell_index) const {
	int x = cell_index / this->get_number_of_micro_cells_high();
	int y = cell_index % this->get_number_of_micro_cells_high();
	return {x, y};
}

bool Param::is_in_bounds(MicroCellPosition position) const {
	return position.x >= 0
		&& position.y >= 0
		&& position.x < this->get_number_of_micro_cells_wide()
		&& position.y < this->get_number_of_micro_cells_high();
}

int Param::get_micro_cell_index_from_position(MicroCellPosition position) const {
	int x = (position.x + this->get_number_of_micro_cells_wide()) % this->get_number_of_micro_cells_wide();
	int y = (position.y + this->get_number_of_micro_cells_high()) % this->get_number_of_micro_cells_high();
	return x * this->get_number_of_micro_cells_high() + y;
}
