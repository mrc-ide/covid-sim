#include "Param.h"

MicroCellPosition Param::get_micro_cell_position_from_cell_index(int cell_index) const {
	int x = cell_index / this->total_microcells_high_;
	int y = cell_index % this->total_microcells_high_;
	return {x, y};
}

bool Param::is_in_bounds(MicroCellPosition const& position) const {
	return position.x >= 0
		&& position.y >= 0
		&& position.x < this->total_microcells_wide_
		&& position.y < this->total_microcells_high_;
}

int Param::get_micro_cell_index_from_position(MicroCellPosition const& position) const {
	int x = (position.x + this->total_microcells_wide_) % this->total_microcells_wide_;
	int y = (position.y + this->total_microcells_high_) % this->total_microcells_high_;
	return x * this->total_microcells_high_ + y;
}
