#pragma once

struct Location {
	float x;
	float y;

	Location();
	Location(float x, float y);
	Location(double x, double y);

	float distance_to(Location other) const;
	float distance_squared_to(Location other) const;
};

