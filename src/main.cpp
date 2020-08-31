#include <iostream>
#include <random>
#include <queue>
#include <unordered_map>
#include <list>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>

#include "extern/poisson_disk_sampling.h"
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
unsigned char REDCOLOR[3] = {255, 0, 0};
unsigned char BLU[3] = {0, 0, 255};
unsigned char DARKGREEN[3] = {0, 100, 0};

void draw_filled_circle(int x0, int y0, int radius, unsigned char *image, int width, int height, int nchannels, unsigned char *color)
{
	int x = radius;
	int y = 0;
	int xchange = 1 - (radius << 1);
	int ychange = 0;
	int err = 0;

	while (x >= y) {
		for (int i = x0 - x; i <= x0 + x; i++) {
			plot(i, y0 + y, image, width, height, nchannels, color);
			plot(i, y0 - y, image, width, height, nchannels, color);
		}
		for (int i = x0 - y; i <= x0 + y; i++) {
			plot(i, y0 + x, image, width, height, nchannels, color);
			plot(i, y0 - x, image, width, height, nchannels, color);
		}

		y++;
		err += ychange;
		ychange += 2;
		if (((err << 1) + xchange) > 0) {
			x--;
			err += xchange;
			xchange += 2;
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

void plant_trees(std::vector<glm::vec2> &tree_points, const Sitemap *map)
{
	float radius = 16.F;
	auto mmin = std::array<float, 2>{{SITE_AREA.min.x, SITE_AREA.min.y}};
	auto mmax = std::array<float, 2>{{SITE_AREA.max.x, SITE_AREA.max.y}};

	std::vector<std::array<float, 2>> candidates = thinks::PoissonDiskSampling(radius, mmin, mmax, 30, map->seed);
	std::vector<glm::vec2> locations;
	for (const auto &point : candidates) {
		locations.push_back(glm::vec2(point[0], point[1]));
	}

	struct byteimage image = blank_byteimage(1, 1024, 1024);
	unsigned char color[] = {255};
	for (const auto &d : map->districts) {
		if (d.radius > 5) {
			glm::vec2 a = {round(d.center.x), round(d.center.y)};
			for (const auto &s : d.sections) {
				glm::vec2 b = {round(s->j0->position.x), round(s->j0->position.y)};
				glm::vec2 c = {round(s->j1->position.x), round(s->j1->position.y)};
				// if street
				if (s->j0->street || s->j1->street) {
					b = a + (0.9F * (b-a));
					c = a + (0.9F * (c-a));
				}
				draw_triangle(a, b, c, image.data, image.width, image.height, image.nchannels, color);
			}
		}
	}

	for (const auto point : locations) {
		if (sample_byteimage(point.x, point.y, RED, &image) > 0.5F) {
			tree_points.push_back(point);
		}
	}

	delete_byteimage(&image);
}

void print_site(const Sitemap *map) 
{
	struct byteimage image = blank_byteimage(3, 1024, 1024);

	std::mt19937 gen(map->seed);
	std::uniform_real_distribution<float> dist(0.f, 1.f);

	unsigned char color[3] = {255, 255, 255};
	for (const auto &d : map->districts) {
		glm::vec2 a = {round(d.center.x), round(d.center.y)};
		glm::vec3 greencolor = {0.78f, 1.f, 0.51f};
		glm::vec3 yellowcolor = {0.91f, 1.f, 0.49f};
		//float gradient = 1.f - (d.radius / 8.f);
		float gradient = 0.75f;
		if (d.radius < 2) {
			gradient = 1.f;
		}
		color[0] = 200 * gradient;
		color[1] = 200 * gradient;
		color[2] = 200 * gradient;
		if (d.radius < 3) {
			color[0] = 181;
			color[1] = 204;
			color[2] = 135;
		}
		if (d.radius > 2 && d.radius < 6) {
			glm::vec3 rgb = glm::mix(yellowcolor, greencolor, dist(gen));
			color[0] = 200*rgb.x;
			color[1] = 200*rgb.y;
			color[2] = 200*rgb.z;
		}
		if (d.radius > 5)  {
			color[0] = 100;
			color[1] = 200;
			color[2] = 0;
		}
		float area = 0;
		for (const auto &s : d.sections) {
			glm::vec2 b = {round(s->j0->position.x), round(s->j0->position.y)};
			glm::vec2 c = {round(s->j1->position.x), round(s->j1->position.y)};
			draw_triangle(a, b, c, image.data, image.width, image.height, image.nchannels, color);
			area += triangle_area(a, b, c);
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

	// draw houses
	//const float AREA_SMALL_HOUSE = 25.F;
	const float SMALL_HOUSE_WIDTH = 5.F;
	const float SMALL_HOUSE_HEIGHT = 10.F;
	const float MEDIUM_HOUSE_WIDTH = 10.F;
	const float MEDIUM_HOUSE_HEIGHT = 20.F;
	const float LARGE_HOUSE_WIDTH = 20.F;
	const float LARGE_HOUSE_HEIGHT = 40.F;
	//const float AREA_LARGE_HOUSE = 200.F;
	for (const auto &parc : map->parcels) {
		float parc_width_front = glm::distance(parc.frontleft, parc.frontright);
		float parc_width_back = glm::distance(parc.backleft, parc.backright);
		float parc_height_left = glm::distance(parc.frontleft, parc.backright);
		float parc_height_right = glm::distance(parc.frontright, parc.backleft);
		float theta = atan2(parc.direction.y, parc.direction.x);
		glm::vec2 perpleft = {-parc.direction.y, parc.direction.x};
		glm::vec2 perpright = {parc.direction.y, -parc.direction.x};
		if ((parc_width_front > LARGE_HOUSE_WIDTH && parc_width_back > LARGE_HOUSE_WIDTH) && (parc_height_left > LARGE_HOUSE_HEIGHT && parc_height_right > LARGE_HOUSE_HEIGHT)) {
			//draw_filled_circle(parc.centroid.x, parc.centroid.y, 12, image.data, image.width, image.height, image.nchannels, PURPLE);
			glm::vec2 a = parc.centroid + (0.5f*LARGE_HOUSE_HEIGHT)*parc.direction + (0.5f*LARGE_HOUSE_WIDTH)*perpleft;
			glm::vec2 b = parc.centroid + (0.5f*LARGE_HOUSE_HEIGHT)*parc.direction + (0.5f*LARGE_HOUSE_WIDTH)*perpright;
			glm::vec2 c = parc.centroid + (0.5f*LARGE_HOUSE_HEIGHT)*(-parc.direction) + (0.5f*LARGE_HOUSE_WIDTH)*perpleft;
			glm::vec2 d = parc.centroid + (0.5f*LARGE_HOUSE_HEIGHT)*(-parc.direction) + (0.5f*LARGE_HOUSE_WIDTH)*perpright;
			draw_thick_line(a.x, a.y, b.x, b.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
			draw_thick_line(a.x, a.y, c.x, c.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
			draw_thick_line(b.x, b.y, d.x, d.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
			draw_thick_line(c.x, c.y, d.x, d.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
			draw_line(parc.centroid.x, parc.centroid.y, parc.centroid.x+0.5f*LARGE_HOUSE_HEIGHT*parc.direction.x, parc.centroid.y+0.5f*LARGE_HOUSE_HEIGHT*parc.direction.y, image.data, image.width, image.height, image.nchannels, REDCOLOR);
			//draw_triangle(parc.centroid, a, b, image.data, image.width, image.height, image.nchannels, REDCOLOR);
		} else if (parc_width_front > MEDIUM_HOUSE_WIDTH && parc_width_back > MEDIUM_HOUSE_WIDTH && parc_height_left > MEDIUM_HOUSE_HEIGHT && parc_height_right > MEDIUM_HOUSE_HEIGHT) {
			glm::vec2 a = parc.centroid + (0.5f*MEDIUM_HOUSE_HEIGHT)*parc.direction + (0.5f*MEDIUM_HOUSE_WIDTH)*perpleft;
			glm::vec2 b = parc.centroid + (0.5f*MEDIUM_HOUSE_HEIGHT)*parc.direction + (0.5f*MEDIUM_HOUSE_WIDTH)*perpright;
			glm::vec2 c = parc.centroid + (0.5f*MEDIUM_HOUSE_HEIGHT)*(-parc.direction) + (0.5f*MEDIUM_HOUSE_WIDTH)*perpleft;
			glm::vec2 d = parc.centroid + (0.5f*MEDIUM_HOUSE_HEIGHT)*(-parc.direction) + (0.5f*MEDIUM_HOUSE_WIDTH)*perpright;
			draw_thick_line(a.x, a.y, b.x, b.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
			draw_thick_line(a.x, a.y, c.x, c.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
			draw_thick_line(b.x, b.y, d.x, d.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
			draw_thick_line(c.x, c.y, d.x, d.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
			draw_line(parc.centroid.x, parc.centroid.y, parc.centroid.x+0.5f*MEDIUM_HOUSE_HEIGHT*parc.direction.x, parc.centroid.y+0.5f*MEDIUM_HOUSE_HEIGHT*parc.direction.y, image.data, image.width, image.height, image.nchannels, REDCOLOR);
		} else if (parc_width_front > SMALL_HOUSE_WIDTH && parc_width_back > SMALL_HOUSE_WIDTH && parc_height_left > SMALL_HOUSE_HEIGHT && parc_height_right > SMALL_HOUSE_HEIGHT) {
			//draw_filled_circle(parc.centroid.x, parc.centroid.y, 3, image.data, image.width, image.height, image.nchannels, PURPLE);
			glm::vec2 a = parc.centroid + (0.5f*SMALL_HOUSE_HEIGHT)*parc.direction + (0.5f*SMALL_HOUSE_WIDTH)*perpleft;
			glm::vec2 b = parc.centroid + (0.5f*SMALL_HOUSE_HEIGHT)*parc.direction + (0.5f*SMALL_HOUSE_WIDTH)*perpright;
			glm::vec2 c = parc.centroid + (0.5f*SMALL_HOUSE_HEIGHT)*(-parc.direction) + (0.5f*SMALL_HOUSE_WIDTH)*perpleft;
			glm::vec2 d = parc.centroid + (0.5f*SMALL_HOUSE_HEIGHT)*(-parc.direction) + (0.5f*SMALL_HOUSE_WIDTH)*perpright;
			draw_thick_line(a.x, a.y, b.x, b.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
			draw_thick_line(a.x, a.y, c.x, c.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
			draw_thick_line(b.x, b.y, d.x, d.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
			draw_thick_line(c.x, c.y, d.x, d.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
			draw_line(parc.centroid.x, parc.centroid.y, parc.centroid.x+0.5f*SMALL_HOUSE_HEIGHT*parc.direction.x, parc.centroid.y+0.5f*SMALL_HOUSE_HEIGHT*parc.direction.y, image.data, image.width, image.height, image.nchannels, REDCOLOR);
		}
	}

	// print farms
	std::bernoulli_distribution bern(0.25f);
	std::bernoulli_distribution tree(0.5f);
	for (const auto &district : map->districts) {
		if (district.radius > 2 && district.radius < 6) {
			if (bern(gen) == true) {
				float max = 0.f;
				struct section *longest = nullptr;
				for (const auto &sect : district.sections) {
					float dist = glm::distance(sect->j0->position, sect->j1->position);
					// also plant trees
					if (tree(gen) == true && dist > 24.F) {
						glm::vec2 start = district.center + 0.9f * (sect->j0->position - district.center);
						glm::vec2 end = district.center + 0.9f * (sect->j1->position - district.center);
						glm::vec2 mid = segment_midpoint(start, end);
						draw_filled_circle(start.x, start.y, 2, image.data, image.width, image.height, image.nchannels, DARKGREEN);
						draw_filled_circle(mid.x, mid.y, 2, image.data, image.width, image.height, image.nchannels, DARKGREEN);
						draw_filled_circle(end.x, end.y, 2, image.data, image.width, image.height, image.nchannels, DARKGREEN);
					}
					if (dist > max) {
						max = dist;
						longest = sect;
					}
				}
				if (longest != nullptr) {
					glm::vec2 dir = glm::normalize(segment_midpoint(longest->j0->position, longest->j1->position) - district.center);
					glm::vec2 perpleft = {-dir.y, dir.x};
					glm::vec2 perpright = {dir.y, -dir.x};

					glm::vec2 a = district.center + (0.5f*MEDIUM_HOUSE_HEIGHT)*dir + (0.5f*MEDIUM_HOUSE_WIDTH)*perpleft;
					glm::vec2 b = district.center + (0.5f*MEDIUM_HOUSE_HEIGHT)*dir + (0.5f*MEDIUM_HOUSE_WIDTH)*perpright;
					glm::vec2 c = district.center + (0.5f*MEDIUM_HOUSE_HEIGHT)*(-dir) + (0.5f*MEDIUM_HOUSE_WIDTH)*perpleft;
					glm::vec2 d = district.center + (0.5f*MEDIUM_HOUSE_HEIGHT)*(-dir) + (0.5f*MEDIUM_HOUSE_WIDTH)*perpright;
					draw_thick_line(a.x, a.y, b.x, b.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
					draw_thick_line(a.x, a.y, c.x, c.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
					draw_thick_line(b.x, b.y, d.x, d.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
					draw_thick_line(c.x, c.y, d.x, d.y, 1, image.data, image.width, image.height, image.nchannels, GRAY);
				}
			}
		}
	}

	// print forests
	std::vector<glm::vec2> tree_points;
	plant_trees(tree_points, map);
	for (const auto &point : tree_points) {
		draw_filled_circle(point.x, point.y, 2, image.data, image.width, image.height, image.nchannels, DARKGREEN);
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
	//long seed = 4207165172587737803; // inwardiness
	//long seed = 6570940619616775359; // infinite

	std::cout << seed << std::endl;

	Sitemap sitemap = {seed, SITE_AREA};

	print_site(&sitemap);

	return 0;
}
