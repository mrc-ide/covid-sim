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