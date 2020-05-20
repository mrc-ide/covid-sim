#include "Param.h"

int Param::get_number_of_micro_cells_wide() const {
	return this->ncw * this->NMCL;
}

int Param::get_number_of_micro_cells_high() const {
	return this->nch * this->NMCL;
}