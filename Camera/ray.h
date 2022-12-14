#pragma once

#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class ray
{
public:
	ray() {}
	ray(const point3& origin, const vec3& direction)
		: orig(origin), dir(direction)
	{}

	point3 origin() const { return orig; } // A
	vec3 direction() const { return dir; } // b

	point3 at(double t) const 
	{
		return orig + t * dir; //A + t * b  
	}
private:
	point3 orig;
	vec3 dir;
};

#endif 