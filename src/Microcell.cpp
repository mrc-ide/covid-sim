#include "Microcell.h"

double Microcell::distance_squared_to(Microcell *other) const {
	double x, y;
	int l, m;

	l = (int)(a - Mcells);
	m = (int)(b - Mcells);
	if (P.DoUTM_coords)
		return dist2UTM(P.in_microcells_.width_ * fabs((double)(l / P.get_number_of_micro_cells_high())), P.in_microcells_.height_ * fabs((double)(l % P.get_number_of_micro_cells_high())),
		                P.in_microcells_.width_ * fabs((double)(m / P.get_number_of_micro_cells_high())), P.in_microcells_.height_ * fabs((double)(m % P.get_number_of_micro_cells_high())));
	else
	{
		x = P.in_microcells_.width_ * fabs((double)(l / P.get_number_of_micro_cells_high() - m / P.get_number_of_micro_cells_high()));
		y = P.in_microcells_.height_ * fabs((double)(l % P.get_number_of_micro_cells_high() - m % P.get_number_of_micro_cells_high()));
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.in_degrees_.width_ * 0.5) x = P.in_degrees_.width_ - x;
			if (y > P.in_degrees_.height_ * 0.5) y = P.in_degrees_.height_ - y;
		}
		return x * x + y * y;
	}
}
