#pragma once

#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

using std::sqrt;

class vec3
{
public:
	vec3() : e{0,0,0} {}
	vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}

	double x() const { return e[0]; }
	double y() const { return e[1]; }
	double z() const { return e[2]; }

	vec3 operator-() const { return vec3(-e[0], -e[1],-e[2]); } // унарный минус
	double operator[](int i) const { return e[i]; }
	double& operator[](int i) { return e[i]; }

	vec3& operator+=(const vec3 &v)
	{
		e[0] += v[0];
		e[1] += v[1];
		e[2] += v[2];
		return *this;
	}

	vec3& operator*=(const double t)
	{
		e[0] *= t;
		e[1] *= t;
		e[2] *= t;
		return *this;
	}

	vec3& operator/=(const double t)
	{
		return *this *= 1 / t; //как это работает??? (мы перегрузили *= и поэтому можем его использовать)
	}

	double length() const
	{
		return sqrt(length_squared());
	}

	double length_squared() const
	{
		return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
	}

	bool near_zero() const {
		// Return true if the vector is close to zero in all dimensions.
		const auto s = 1e-8;
		return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
	}

	inline friend std::ostream& operator<<(std::ostream& out, const vec3& v)
	{
		return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
	}

	inline friend vec3 operator+(const vec3& u, const vec3& v) //а почему не vec& возращаем
	{
		return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
	}

	inline friend vec3 operator-(const vec3& u, const vec3& v)
	{
		return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
	}

	inline friend vec3 operator*(const vec3& u, const vec3& v)
	{
		return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
	}

	inline friend vec3 operator*(const double t, const vec3& v)
	{
		return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
	}

	inline friend vec3 operator*(const vec3& v, const double t)
	{
		return t * v;
	}

	inline friend vec3 operator/(vec3 v, const double t) // почему не vec& v
	{
		return (1/t) * v;
	}

	static inline double dot(const vec3& u, const vec3& v)
	{
		return u.e[0] * v.e[0]
			+ u.e[1] * v.e[1]
			+ u.e[2] * v.e[2];
	}

	static inline vec3 cross(const vec3& u, const vec3& v)
	{
		return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
					u.e[2] * v.e[0] - u.e[0] * v.e[2],
					u.e[0] * v.e[1] - u.e[1] * v.e[0]);
	}



	

private:
	double e[3];
};

inline vec3 unit_vector(vec3 v) //единичний (направляющий) вектор
{
	return v / v.length();
}


inline vec3 random()
{
	return vec3(random_double(), random_double(), random_double());
}

inline vec3 random(double min, double max)
{
	return vec3(random_double(min, max), random_double(min, max), random_double(min, max));
}

inline vec3 random_in_unit_sphere()
{
	while (true)
	{
		auto p = random(-1, 1);
		if (p.length_squared() >= 1) continue;
		return p;
	}
}

inline vec3 random_unit_vector()
{
	return unit_vector(random_in_unit_sphere());
}

vec3 random_in_hemisphere(const vec3& normal)
{
	vec3 in_unit_sphere = random_in_unit_sphere();
	if (vec3::dot(in_unit_sphere, normal) > 0.0)
	{
		return in_unit_sphere;
	}
	else 
	{
		return -in_unit_sphere;
	}
}

vec3 random_in_unit_disk()
{
	while (true)
	{
		auto p = vec3(random_double(-1, 1), random_double(-1, 1), 0);
		if (p.length_squared() >= 1) continue;
		return p;
	}
}

vec3 reflect(const vec3& v, const vec3& n) {
	return v - 2 * vec3::dot(v, n) * n;
}

vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
	auto cos_theta = fmin(vec3::dot(-uv, n), 1.0);
	vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
	vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
	return (r_out_perp + r_out_parallel);
}
// Type aliases for vec3 
using point3 = vec3;
using color = vec3;

#endif
