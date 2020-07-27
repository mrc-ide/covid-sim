#include <cstdlib>
#include <cmath>

#include "Constants.h"
#include "Dist.h"
#include "Param.h"

#include "Model.h"

double sinx[DEGREES_PER_TURN + 1], cosx[DEGREES_PER_TURN + 1], asin2sqx[1001];

//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****
//// **** DISTANCE FUNCTIONS (return distance-squared, which is input for every Kernel function)
//// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// **** //// ****

double periodic_xy(double x, double y) {
	if (P.DoPeriodicBoundaries)
	{
		if (x > P.in_degrees_.width * 0.5) x = P.in_degrees_.width - x;
		if (y > P.in_degrees_.height * 0.5) y = P.in_degrees_.height - y;
	}
	return x * x + y * y;
}

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
	yt = fabs(y1 + P.SpatialBoundingBox.bottom_left().y);
	yi = floor(yt);
	cy1 = yt - yi;
	cy1 = (1 - cy1) * cosx[((int)yi)] + cy1 * cosx[((int)yi) + 1];
	yt = fabs(y2 + P.SpatialBoundingBox.bottom_left().y);
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
	double x, y;

	if (P.DoUTM_coords)
		return dist2UTM(Households[a->hh].loc.x, Households[a->hh].loc.y, Households[b->hh].loc.x, Households[b->hh].loc.y);
	else
	{
		x = fabs(Households[a->hh].loc.x - Households[b->hh].loc.x);
		y = fabs(Households[a->hh].loc.y - Households[b->hh].loc.y);
		return periodic_xy(x, y);
	}
}
double dist2_cc(Cell* a, Cell* b)
{
	double x, y;
	int l, m;

	l = (int)(a - Cells);
	m = (int)(b - Cells);
	if (P.DoUTM_coords)
		return dist2UTM(P.in_cells_.width * fabs((double)(l / P.nch)), P.in_cells_.height * fabs((double)(l % P.nch)),
			P.in_cells_.width * fabs((double)(m / P.nch)), P.in_cells_.height * fabs((double)(m % P.nch)));
	else
	{
		x = P.in_cells_.width * fabs((double)(l / P.nch - m / P.nch));
		y = P.in_cells_.height * fabs((double)(l % P.nch - m % P.nch));
		return periodic_xy(x, y);
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
		if (P.in_cells_.width * ((double)abs(m / P.nch - l / P.nch)) > PI)
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
		return dist2UTM(P.in_cells_.width * fabs((double)(i / P.nch)), P.in_cells_.height * fabs((double)(i % P.nch)),
			P.in_cells_.width * fabs((double)(j / P.nch)), P.in_cells_.height * fabs((double)(j % P.nch)));
	}
	else
	{
		if ((P.DoPeriodicBoundaries) && (P.in_cells_.width * ((double)abs(m / P.nch - l / P.nch)) > P.in_degrees_.width * 0.5))
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
		if ((P.DoPeriodicBoundaries) && (P.in_degrees_.height * ((double)abs(m % P.nch - l % P.nch)) > P.in_degrees_.height * 0.5))
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
		x = P.in_cells_.width * fabs((double)(i / P.nch - j / P.nch));
		y = P.in_cells_.height * fabs((double)(i % P.nch - j % P.nch));
		return periodic_xy(x, y);
	}
}
double dist2_mm(Microcell* a, Microcell* b)
{
	double x, y;
	int l, m;

	l = (int)(a - Mcells);
	m = (int)(b - Mcells);
	if (P.DoUTM_coords)
		return dist2UTM(P.in_microcells_.width * fabs((double)(l / P.total_microcells_high_)), P.in_microcells_.height * fabs((double)(l % P.total_microcells_high_)),
			P.in_microcells_.width * fabs((double)(m / P.total_microcells_high_)), P.in_microcells_.height * fabs((double)(m % P.total_microcells_high_)));
	else
	{
		x = P.in_microcells_.width * fabs((double)(l / P.total_microcells_high_ - m / P.total_microcells_high_));
		y = P.in_microcells_.height * fabs((double)(l % P.total_microcells_high_ - m % P.total_microcells_high_));
		return periodic_xy(x, y);
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
		return periodic_xy(x, y);
	}
}
