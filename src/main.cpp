#include <iostream>
#include <random>
#include <queue>
#include <unordered_map>
#include <list>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "extern/stb_image_write.h"

#include "geom.h"
#include "imp.h"
#include "voronoi.h"
#include "sitemap.h"

const struct rectangle SITE_AREA = { .min = {0.F, 0.F}, .max = {1024.F, 1024.F} };
const glm::vec2 CENTER = {512.F, 512.F};
unsigned char PURPLE[3] = {255, 0, 255};
unsigned char ORANGE[3] = {255, 150, 0};
unsigned char GRAY[3] = {112, 112, 112};
unsigned char LIGHTGRAY[3] = {150, 150, 150};
unsigned char BLACK[3] = {0, 0, 0};
unsigned char GRN[3] = {0, 255, 0};
unsigned char YELLOW[3] = {255, 255, 0};

void draw_filled_circle(int x0, int y0, int radius, unsigned char *image, int width, int height, int nchannels, unsigned char *color)
{
	int x = radius;
	int y = 0;
	int xChange = 1 - (radius << 1);
	int yChange = 0;
	int radiusError = 0;

	while (x >= y)
	{
	for (int i = x0 - x; i <= x0 + x; i++)
	{
	plot(i, y0 + y, image, width, height, nchannels, color);
	plot(i, y0 - y, image, width, height, nchannels, color);
	}
	for (int i = x0 - y; i <= x0 + y; i++)
	{
	plot(i, y0 + x, image, width, height, nchannels, color);
	plot(i, y0 - x, image, width, height, nchannels, color);
	}

	y++;
	radiusError += yChange;
	yChange += 2;
	if (((radiusError << 1) + xChange) > 0)
	{
	x--;
	radiusError += xChange;
	xChange += 2;
	}
	}
}

void draw_thick_line(int x0, int y0, int x1, int y1, int radius, unsigned char *image, int width, int height, int nchannels, unsigned char *color)
{
	int dx =  abs(x1-x0), sx = x0<x1 ? 1 : -1;
	int dy = -abs(y1-y0), sy = y0<y1 ? 1 : -1;
	int err = dx+dy, e2; // error value e_xy

	for(;;) {
		draw_filled_circle(x0,y0, radius, image, width, height, nchannels, color);
		if (x0==x1 && y0==y1) { break; }
		e2 = 2*err;
		if (e2 >= dy) { err += dy; x0 += sx; } // e_xy+e_x > 0
		if (e2 <= dx) { err += dx; y0 += sy; } // e_xy+e_y < 0
	}
}

void print_site(const Sitemap *map) 
{
	struct byteimage image = blank_byteimage(3, 1024, 1024);

	unsigned char color[3] = {255, 255, 255};
	for (const auto &d : map->districts) {
		glm::vec2 a = {round(d.center.x), round(d.center.y)};
		//float gradient = 1.f - (d.radius / 8.f);
		float gradient = 0.75f;
		if (d.radius < 2) {
			gradient = 1.f;
		}
		color[0] = 200 * gradient;
		color[1] = 200 * gradient;
		color[2] = 200 * gradient;
		for (const auto &s : d.sections) {
			glm::vec2 b = {round(s->j0->position.x), round(s->j0->position.y)};
			glm::vec2 c = {round(s->j1->position.x), round(s->j1->position.y)};
			draw_triangle(a, b, c, image.data, image.width, image.height, image.nchannels, color);
		}
		//draw_filled_circle(a.x, a.y, 1, image.data, image.width, image.height, image.nchannels, BLACK);
	}
	/*
	for (const auto &sect : map->sections) {
		draw_thick_line(sect.j0->position.x, sect.j0->position.y, sect.j1->position.x, sect.j1->position.y, 1, image.data, image.width, image.height, image.nchannels, GRN);
	}
	*/

	for (const auto &sect : map->sections) {
		if (sect.j0->radius < 3 && sect.j1->radius < 3) {
			if (sect.j0->border == false && sect.j1->border == false) {
				draw_thick_line(sect.j0->position.x, sect.j0->position.y, sect.j1->position.x, sect.j1->position.y, 3, image.data, image.width, image.height, image.nchannels, ORANGE);
			}
		}
	}

	 for (const auto &d : map->districts) {
		if (d.radius < 2) {
			float max = 0.f;
			const struct section *longest;
			for (const auto &s : d.sections) {
				if (s->j0-> radius < 7 || s->j1->radius < 7) {
					glm::vec2 mid = segment_midpoint(s->j0->position, s->j1->position);
					float dist = glm::distance(s->j0->position, s->j1->position);
					if (dist > max) {
						max = dist;
						longest = s;
					}
				}
			}
			glm::vec2 mad = segment_midpoint(longest->j0->position, longest->j1->position);
			glm::vec2 tjunction = glm::normalize(d.center - mad);
			draw_thick_line(d.center.x, d.center.y, mad.x, mad.y, 3, image.data, image.width, image.height, image.nchannels, ORANGE);
			const struct section *coolest;
			float mindot = 1.f;
			for (const auto &s : d.sections) {
				if (s->j0-> radius < 8 || s->j1->radius < 8) {
					glm::vec2 mid = segment_midpoint(s->j0->position, s->j1->position);
					glm::vec2 cool = glm::normalize(d.center - mid);
					float dot = glm::dot(tjunction, cool);
					if (dot < mindot) {
						mindot = dot;
						coolest = s;
					}
				}
			}
			glm::vec2 mod = segment_midpoint(coolest->j0->position, coolest->j1->position);
			draw_thick_line(d.center.x, d.center.y, mod.x, mod.y, 3, image.data, image.width, image.height, image.nchannels, ORANGE);
		}
	}

	for (const auto &way : map->highways) {
		draw_thick_line(way.P0.x, way.P0.y, way.P1.x, way.P1.y, 4, image.data, image.width, image.height, image.nchannels, ORANGE);
	}


	for (const auto &sect : map->sections) {
		if (sect.wall) {
			draw_thick_line(sect.d0->center.x, sect.d0->center.y, sect.d1->center.x, sect.d1->center.y, 6, image.data, image.width, image.height, image.nchannels, GRAY);
		}
	}
	for (const auto &d : map->districts) {
		if (d.wall) {
			draw_filled_circle(d.center.x, d.center.y, 12, image.data, image.width, image.height, image.nchannels, GRAY);
		}
	}


	draw_filled_circle(map->core->center.x, map->core->center.y, 4, image.data, image.width, image.height, image.nchannels, BLACK);

	for (const auto &sect : map->sections) {
		if (sect.gateway) {
			draw_thick_line(sect.d0->center.x, sect.d0->center.y, sect.d1->center.x, sect.d1->center.y, 6, image.data, image.width, image.height, image.nchannels, GRAY);
			glm::vec2 outward = segment_midpoint(sect.d0->center, sect.d1->center);
			glm::vec2 right = sect.d0->center - outward;
			glm::vec2 a = sect.j0->position - right;;
			glm::vec2 b = sect.j0->position + right;;
			glm::vec2 c = sect.j1->position - right;;
			glm::vec2 d = sect.j1->position + right;;
			draw_thick_line(a.x, a.y, b.x, b.y, 1, image.data, image.width, image.height, image.nchannels, YELLOW);
			draw_thick_line(c.x, c.y, d.x, d.y, 1, image.data, image.width, image.height, image.nchannels, YELLOW);
			draw_thick_line(a.x, a.y, c.x, c.y, 1, image.data, image.width, image.height, image.nchannels, YELLOW);
			draw_thick_line(b.x, b.y, d.x, d.y, 1, image.data, image.width, image.height, image.nchannels, YELLOW);
		}
	}

	stbi_flip_vertically_on_write(true);
	stbi_write_png("diagram.png", image.width, image.height, image.nchannels, image.data, image.width*image.nchannels);

	delete_byteimage(&image);
}

int main(int argc, char *argv[])
{
	using namespace glm;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<long> dist;
	long seed = dist(gen);

	std::cout << seed << std::endl;

	Sitemap sitemap = {seed, SITE_AREA};

	print_site(&sitemap);

	return 0;
}
