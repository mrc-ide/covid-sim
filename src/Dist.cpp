#include <stdlib.h>
#include <math.h>

#include "Constants.h"
#include "Dist.h"
#include "Param.h"

double sinx[361], cosx[361], asin2sqx[1001];

//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
//// **** DISTANCE FUNCTIONS (return distance-squared, which is input for every Kernel function)
//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

double dist2UTM(double x1, double y1, double x2, double y2)
{
	double x, y, cy1, cy2, yt, xi, yi;

	x = fabs(x1 - x2) / 2;
	y = fabs(y1 - y2) / 2;
	xi = floor(x);
	yi = floor(y);
	x -= xi;
	y -= yi;
	x = (1 - x) * sinx[(int)xi] + x * sinx[((int)xi) + 1];
	y = (1 - y) * sinx[(int)yi] + y * sinx[((int)yi) + 1];
	yt = fabs(y1 + P.SpatialBoundingBox[1]);
	yi = floor(yt);
	cy1 = yt - yi;
	cy1 = (1 - cy1) * cosx[((int)yi)] + cy1 * cosx[((int)yi) + 1];
	yt = fabs(y2 + P.SpatialBoundingBox[1]);
	yi = floor(yt);
	cy2 = yt - yi;
	cy2 = (1 - cy2) * cosx[((int)yi)] + cy2 * cosx[((int)yi) + 1];
	x = fabs(1000 * (y * y + x * x * cy1 * cy2));
	xi = floor(x);
	x -= xi;
	y = (1 - x) * asin2sqx[((int)xi)] + x * asin2sqx[((int)xi) + 1];
	return 4 * EARTHRADIUS * EARTHRADIUS * y;
}
double dist2(Person* a, Person* b)
{
	if (P.DoUTM_coords)
		return Households[a->hh].distance_squared_to(Households[b->hh]);
	else
	{
		Vector2<float> delta = (Households[a->hh].loc - Households[b->hh].loc).abs();
		if (P.DoPeriodicBoundaries)
		{
			if (delta.x > P.in_degrees_.width * 0.5) delta.x = P.in_degrees_.width - delta.x;
			if (delta.y > P.in_degrees_.height * 0.5) delta.y = P.in_degrees_.height - delta.y;
		}
		return delta.length_squared();
	}
}

double dist2_cc(Cell* a, Cell* b)
{
	int aIndex = (int)(a - Cells);
	int bIndex = (int)(b - Cells);
	if (P.DoUTM_coords) {
		return dist2UTM(P.in_cells_.width_ * fabs((double) (aIndex / P.nch)),
		                P.in_cells_.height_ * fabs((double) (aIndex % P.nch)),
		                P.in_cells_.width_ * fabs((double) (bIndex / P.nch)),
		                P.in_cells_.height_ * fabs((double) (bIndex % P.nch)));
	}
	else
	{
		double x = P.in_cells_.width_ * fabs((double)(aIndex / P.nch - bIndex / P.nch));
		double y = P.in_cells_.height_ * fabs((double)(aIndex % P.nch - bIndex % P.nch));
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.in_degrees_.width * 0.5) x = P.in_degrees_.width - x;
			if (y > P.in_degrees_.height * 0.5) y = P.in_degrees_.height - y;
		}
		return x * x + y * y;
	}
}
double dist2_cc_min(Cell* a, Cell* b)
{
	double x, y;
	int l, m, i, j;

	l = (int)(a - Cells);
	m = (int)(b - Cells);
	i = l; j = m;
	if (P.DoUTM_coords)
	{
		if (P.in_cells_.width_ * ((double)abs(m / P.nch - l / P.nch)) > PI)
		{
			if (m / P.nch > l / P.nch)
				j += P.nch;
			else if (m / P.nch < l / P.nch)
				i += P.nch;
		}
		else
		{
			if (m / P.nch > l / P.nch)
				i += P.nch;
			else if (m / P.nch < l / P.nch)
				j += P.nch;
		}
		if (m % P.nch > l % P.nch)
			i++;
		else if (m % P.nch < l % P.nch)
			j++;
		return dist2UTM(P.in_cells_.width_ * fabs((double)(i / P.nch)), P.in_cells_.height_ * fabs((double)(i % P.nch)),
			P.in_cells_.width_ * fabs((double)(j / P.nch)), P.in_cells_.height_ * fabs((double)(j % P.nch)));
	}
	else
	{
		if ((P.DoPeriodicBoundaries) && (P.in_cells_.width_ * ((double)abs(m / P.nch - l / P.nch)) > P.in_degrees_.width_ * 0.5))
		{
			if (m / P.nch > l / P.nch)
				j += P.nch;
			else if (m / P.nch < l / P.nch)
				i += P.nch;
		}
		else
		{
			if (m / P.nch > l / P.nch)
				i += P.nch;
			else if (m / P.nch < l / P.nch)
				j += P.nch;
		}
		if ((P.DoPeriodicBoundaries) && (P.in_degrees_.height_ * ((double)abs(m % P.nch - l % P.nch)) > P.in_degrees_.height_ * 0.5))
		{
			if (m % P.nch > l % P.nch)
				j++;
			else if (m % P.nch < l % P.nch)
				i++;
		}
		else
		{
			if (m % P.nch > l % P.nch)
				i++;
			else if (m % P.nch < l % P.nch)
				j++;
		}
		x = P.in_cells_.width_ * fabs((double)(i / P.nch - j / P.nch));
		y = P.in_cells_.height_ * fabs((double)(i % P.nch - j % P.nch));
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.in_degrees_.width_ * 0.5) x = P.in_degrees_.width_ - x;
			if (y > P.in_degrees_.height_ * 0.5) y = P.in_degrees_.height_ - y;
		}
		return x * x + y * y;
	}
}
double dist2_mm(Microcell* a, Microcell* b)
{
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

double dist2_raw(double ax, double ay, double bx, double by)
{
	double x, y;

	if (P.DoUTM_coords)
		return dist2UTM(ax, ay, bx, by);
	else
	{
		x = fabs(ax - bx);
		y = fabs(ay - by);
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.in_degrees_.width_ * 0.5) x = P.in_degrees_.width_ - x;
			if (y > P.in_degrees_.height_ * 0.5) y = P.in_degrees_.height_ - y;
		}
		return x * x + y * y;
	}
}
