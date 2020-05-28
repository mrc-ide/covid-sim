#pragma once

#include "Country.h"
#include "Geometry/Vector2.h"

namespace Models {
/**
 * @brief Represents an institution that people may belong to.
 *
 * PLACE be an elementary school, high schools, universities, workplaces etc. Places
 * belong to a microcell (and therefore have a spatial location). Places may have state
 * (i.e., closed or open). Mechanisms exist for absenteeism tracking (but are not currently used).
 * The `members` array lists all individuals who belong to a place.
 * Places can have different groups (to model differential interaction strengths between groups
 * in the same place).
 */
	struct Place {
		int n, mcell;
		unsigned short int ng, treat, control_trig, country;
		unsigned short int close_start_time, close_end_time, treat_end_time;
		unsigned short int *AvailByAge;
		unsigned short int Absent[MAX_ABSENT_TIME], AbsentLastUpdateTime;
		Geometry::Vector2<float> loc;
		float ProbClose;
		int *group_start, *group_size, *members;
	};
}