#include <cmath>

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

	x = std::fabs(x1 - x2) / 2;
	y = std::fabs(y1 - y2) / 2;
	xi = floor(x);
	yi = floor(y);
	x -= xi;
	y -= yi;
	x = (1 - x) * sinx[(int)xi] + x * sinx[((int)xi) + 1];
	y = (1 - y) * sinx[(int)yi] + y * sinx[((int)yi) + 1];
	yt = std::fabs(y1 + P.SpatialBoundingBox[1]);
	yi = floor(yt);
	cy1 = yt - yi;
	cy1 = (1 - cy1) * cosx[((int)yi)] + cy1 * cosx[((int)yi) + 1];
	yt = std::fabs(y2 + P.SpatialBoundingBox[1]);
	yi = floor(yt);
	cy2 = yt - yi;
	cy2 = (1 - cy2) * cosx[((int)yi)] + cy2 * cosx[((int)yi) + 1];
	x = std::fabs(1000 * (y * y + x * x * cy1 * cy2));
	xi = floor(x);
	x -= xi;
	y = (1 - x) * asin2sqx[((int)xi)] + x * asin2sqx[((int)xi) + 1];
	return 4 * EARTHRADIUS * EARTHRADIUS * y;
}
double dist2(person* a, person* b)
{
	double x, y;

	if (P.DoUTM_coords)
		return dist2UTM(Households[a->hh].loc_x, Households[a->hh].loc_y, Households[b->hh].loc_x, Households[b->hh].loc_y);
	else
	{
		x = std::fabs(Households[a->hh].loc_x - Households[b->hh].loc_x);
		y = std::fabs(Households[a->hh].loc_y - Households[b->hh].loc_y);
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.width * 0.5) x = P.width - x;
			if (y > P.height * 0.5) y = P.height - y;
		}
		return x * x + y * y;
	}
}
double dist2_cc(cell* a, cell* b)
{
	double x, y;
	int l, m;

	l = (int)(a - Cells);
	m = (int)(b - Cells);
	if (P.DoUTM_coords)
		return dist2UTM(P.cwidth * std::fabs((double)(l / P.nch)), P.cheight * std::fabs((double)(l % P.nch)),
			P.cwidth * std::fabs((double)(m / P.nch)), P.cheight * std::fabs((double)(m % P.nch)));
	else
	{
		x = P.cwidth * std::fabs((double)(l / P.nch - m / P.nch));
		y = P.cheight * std::fabs((double)(l % P.nch - m % P.nch));
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.width * 0.5) x = P.width - x;
			if (y > P.height * 0.5) y = P.height - y;
		}
		return x * x + y * y;
	}
}
double dist2_cc_min(cell* a, cell* b)
{
	double x, y;
	int l, m, i, j;

	l = (int)(a - Cells);
	m = (int)(b - Cells);
	i = l; j = m;
	if (P.DoUTM_coords)
	{
		if (P.cwidth * std::fabs(m / P.nch - l / P.nch) > PI)
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
		return dist2UTM(P.cwidth * std::fabs((double)(i / P.nch)), P.cheight * std::fabs((double)(i % P.nch)),
			P.cwidth * std::fabs((double)(j / P.nch)), P.cheight * std::fabs((double)(j % P.nch)));
	}
	else
	{
		if ((P.DoPeriodicBoundaries) && (P.cwidth * std::fabs(m / P.nch - l / P.nch) > P.width * 0.5))
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
		if ((P.DoPeriodicBoundaries) && (P.height * std::fabs(m % P.nch - l % P.nch) > P.height * 0.5))
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
		x = P.cwidth * std::fabs((double)(i / P.nch - j / P.nch));
		y = P.cheight * std::fabs((double)(i % P.nch - j % P.nch));
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.width * 0.5) x = P.width - x;
			if (y > P.height * 0.5) y = P.height - y;
		}
		return x * x + y * y;
	}
}
double dist2_mm(microcell* a, microcell* b)
{
	double x, y;
	int l, m;

	l = (int)(a - Mcells);
	m = (int)(b - Mcells);
	if (P.DoUTM_coords)
		return dist2UTM(P.mcwidth * std::fabs((double)(l / P.nmch)), P.mcheight * std::fabs((double)(l % P.nmch)),
			P.mcwidth * std::fabs((double)(m / P.nmch)), P.mcheight * std::fabs((double)(m % P.nmch)));
	else
	{
		x = P.mcwidth * std::fabs((double)(l / P.nmch - m / P.nmch));
		y = P.mcheight * std::fabs((double)(l % P.nmch - m % P.nmch));
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.width * 0.5) x = P.width - x;
			if (y > P.height * 0.5) y = P.height - y;
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
		x = std::fabs(ax - bx);
		y = std::fabs(ay - by);
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.width * 0.5) x = P.width - x;
			if (y > P.height * 0.5) y = P.height - y;
		}
		return x * x + y * y;
	}
}
