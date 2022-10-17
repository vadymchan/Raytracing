#pragma once

#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"

class sphere : public hittable
{
public:
	sphere() {}
	sphere(point3 cen, double r) : center(cen), radius(r) {}

	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec)const override;
private:
	point3 center;
	double radius;
};

bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec)const
{
	vec3 oc = r.origin() - center; //???
	auto a = r.direction().length_squared();
	auto half_b = vec3::dot(oc, r.direction());
	auto c = oc.length_squared() - radius * radius;

	auto discriminant = half_b * half_b - a * c;
	if (discriminant < 0) return false;
	auto sqrtd = sqrt(discriminant);

	//find the nearest root that lies in the aceptable range.
	auto root = (-half_b - sqrtd) / a;
	if (root < t_min || t_max < root) //вышли за пределы
	{
		root = (-half_b + sqrtd) / a; // 2 корень
		if (root < t_min || t_max < root)
		{
			return false;
		}
	}


	//записываем данные где пересекли сферу
	rec.t = root; //???
	rec.p = r.at(rec.t); 
	rec.normal = (rec.p - center) / radius;

	return true;
	
}


#endif