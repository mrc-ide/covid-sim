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
	Vector2<double> pos = Vector2<double>(x1 - x2, y1 - y2).abs() / 2;
	Vector2<double> i = pos.floor();
	pos -= i;
	pos.x = (1 - pos.x) * sinx[(int)i.x] + pos.x * sinx[((int)i.x) + 1];
	pos.y = (1 - pos.y) * sinx[(int)i.y] + pos.y * sinx[((int)i.y) + 1];

	double yt = fabs(y1 + P.SpatialBoundingBox.start.y);
	i.y = floor(yt);

	double cy1 = yt - i.y;
	cy1 = (1 - cy1) * cosx[((int)i.y)] + cy1 * cosx[((int)i.y) + 1];
	yt = fabs(y2 + P.SpatialBoundingBox.start.y);
	i.y = floor(yt);

	double cy2 = yt - i.y;
	cy2 = (1 - cy2) * cosx[((int)i.y)] + cy2 * cosx[((int)i.y) + 1];

	pos.x = fabs(1000 * (pos.length_squared() * cy1 * cy2));
	i.x = floor(pos.x);
	pos.x -= i.x;
	pos.y = (1 - pos.x) * asin2sqx[((int)i.x)] + pos.x * asin2sqx[((int)i.x) + 1];
	return 4 * EARTHRADIUS * EARTHRADIUS * pos.y;
}

double dist2_raw(double ax, double ay, double bx, double by)
{
	if (P.DoUTM_coords)
		return dist2UTM(ax, ay, bx, by);
	else
	{
		double x = fabs(ax - bx);
		double y = fabs(ay - by);
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.in_degrees_.width * 0.5) x = P.in_degrees_.width - x;
			if (y > P.in_degrees_.height * 0.5) y = P.in_degrees_.height - y;
		}
		return x * x + y * y;
	}
}
