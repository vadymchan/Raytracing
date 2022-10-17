#include "rtweekend.h" 

#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"

#include <iostream>

color ray_color(const ray& r, const hittable& world, int depth)
{
	hit_record rec;

	if (depth <= 0)
		return color(0, 0, 0);

	if (world.hit(r, 0.001, infinity, rec)) // нужно игнорировать попадания близкие к 0
	{
		ray scattered;
		color attenuation;
		if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
			return attenuation * ray_color(scattered, world, depth - 1);
		return color(0, 0, 0);
	}
	vec3 unit_direction = unit_vector(r.direction());
	auto t = 0.5 * (unit_direction.y() + 1.0);
	return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}


double hit_sphere(const point3& center, double radius, const ray& r)
{
	vec3 oc = r.origin() - center;					  // A - C
	/*
	
	auto a = vec3::dot(r.direction(), r.direction()); // b * b
	auto b = 2.0 * vec3::dot(oc, r.direction());	  // 2 * (A - C) * b 
	auto c = vec3::dot(oc, oc) - radius * radius;	  // (A - C) * (A - C) - r * r
	auto discriminant = b * b - 4 * a * c;

	*/

	auto a = r.direction().length_squared();
	auto half_b = vec3::dot(oc, r.direction());
	auto c = oc.length_squared() - radius * radius;
	auto discriminant = half_b * half_b - a * c;

	if (discriminant < 0)
	{
		return -1.0;
	}
	else
	{
		//return ((-b - sqrt(discriminant)) / (2.0 * a));
		return (( - half_b - sqrt(discriminant)) / a);
	}
}

//color ray_color(const ray& r)
//{
//	auto t = hit_sphere(point3(0, 0, -1), 0.5, r);
//	if (t > 0.0)
//	{
//		vec3 N = vec3::unit_vector(r.at(t) - vec3(0, 0, -1));
//		return 0.5 * color(N.x() + 1, N.y() + 1, N.z() + 1);
//	}
//	vec3 unit_direction = vec3::unit_vector(r.direction());
//	t = 0.5 * (unit_direction.y() + 1.0); //направление y ( -1.0 < y < 1.0) кстати, подумай над направляющим вектором
//	return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
//}


hittable_list random_scene()
{
	hittable_list world;

	auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
	world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, ground_material));

	for (int a = -11; a < 11; ++a)
	{
		for (int b = -11; b < 11; b++)
		{
			auto choose_mat = random_double();
			point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

			if ((center - point3(4, 0.2, 0)).length() > 0.9)
			{
				shared_ptr<material> sphere_material;

				if (choose_mat < 0.8)
				{
					//diffuse
					auto albedo = random() * random();
					sphere_material = make_shared<lambertian>(albedo);
				}
				else if (choose_mat < 0.95)
				{
					//metal
					auto albedo = random(0.5, 1);
					auto fuzz = random_double(0, 0.5);
					sphere_material = make_shared<metal>(albedo, fuzz);
				}
				else
				{
					//glass
					sphere_material = make_shared<dielectric>(1.5);
				}
				world.add(make_shared<sphere>(center, 0.2, sphere_material));
			}
		}
	}

	auto material1 = make_shared<dielectric>(1.5);
	world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

	auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
	world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

	auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
	world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

	return world;
}

int main()
{

	//Image 

	const auto aspect_ratio = 3.0 / 2.0;
	const int image_width = 1200;
	const int image_height = static_cast<int>(image_width / aspect_ratio);
	const int samples_per_pixel = 500;
	const int max_depth = 50;

	//World

	//auto R = cos(pi / 4); //???
	hittable_list world = random_scene();
	
	/*auto material_left = make_shared<lambertian>(color(0, 0, 1));
	auto material_right = make_shared<lambertian>(color(1, 0, 0));
	
	world.add(make_shared<sphere>(point3(-R, 0, -1), R, material_left));
	world.add(make_shared<sphere>(point3(R, 0, -1), R, material_right));*/

	//auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
	//auto material_center = make_shared<dielectric>(1.5);
	//auto material_left = make_shared<dielectric>(1.5);
	//auto material_right = make_shared<metal>(color(0.8, 0.6, 0.2), 1.0);

	//auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
	//auto material_center = make_shared<lambertian>(color(0.1, 0.2, 0.5));
	//auto material_left = make_shared<dielectric>(1.5);
	//auto material_right = make_shared<metal>(color(0.8, 0.6, 0.2), 1.0);

	auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
	auto material_center = make_shared<lambertian>(color(0.1, 0.2, 0.5));
	auto material_left = make_shared<dielectric>(1.5);
	auto material_right = make_shared<metal>(color(0.8, 0.6, 0.2), 0.0);

	world.add(make_shared<sphere>(point3(0.0, -100.5, -1.0), 100.0, material_ground));
	world.add(make_shared<sphere>(point3(0.0, 0.0, -1.0), 0.5, material_center));
	world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), 0.5, material_left));
	world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), -0.45, material_left)); // это для hollow sphere
	world.add(make_shared<sphere>(point3(1.0, 0.0, -1.0), 0.5, material_right));

	//Camera

	point3 lookfrom(13, 2, 3);
	point3 lookat(0, 0, 0);
	vec3 vup(0, 1, 0);
	//auto dist_to_focus = (lookfrom - lookat).length();
	auto dist_to_focus = 10;
	auto aperture = 0.1;

	camera cam(lookfrom, lookat, vup, 20.0, aspect_ratio, aperture, dist_to_focus);

	//Render

	std::cout << "P3\n" << image_width << ' ' << image_height << "\n256\n";

	for (int j = image_height - 1; j >= 0; --j)
	{
		std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
		for (int i = 0; i < image_width; ++i)
		{
			//color pixel_color(double(i) / (image_width - 1), double(j) / (image_height - 1), 0.5); // r g b
			color pixel_color(0, 0, 0);
			for (int s = 0; s < samples_per_pixel; ++s)
			{
				auto u = double(i + random_double()) / (image_width - 1);
				auto v = double(j + random_double()) / (image_height - 1);
				//ray r(origin, lower_left_corner + u * horizontal + v * vertical - origin);
				ray r = cam.get_ray(u, v);
				pixel_color += ray_color(r, world, max_depth);
			}
			write_color(std::cout, pixel_color, samples_per_pixel);
		}
	}

	std::cerr << "\nDone.\n";

}