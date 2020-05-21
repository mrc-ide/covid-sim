#pragma once

struct Location {
	const float x;
	const float y;

	float distance_to(Location other) const;
	float distance_squared_to(Location other) const;
};

