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

void print_site(const Sitemap *map) 
{
	struct byteimage image = blank_byteimage(3, 1024, 1024);

	unsigned char color[3] = {255, 255, 255};
	for (const auto &d : map->districts) {
		glm::vec2 a = {round(d.center.x), round(d.center.y)};
		float gradient = 1.f - (d.radius / 8.f);
		color[0] = 255 * gradient;
		color[1] = 255 * gradient;
		color[2] = 255 * gradient;
		for (const auto &s : d.sections) {
			glm::vec2 b = {round(s->j0->position.x), round(s->j0->position.y)};
			glm::vec2 c = {round(s->j1->position.x), round(s->j1->position.y)};
			draw_triangle(a, b, c, image.data, image.width, image.height, image.nchannels, color);
		}
	}

	for (const auto &j : map->junctions) {
		if (j.radius == 2) {
			plot(j.position.x, j.position.y, image.data, image.width, image.height, image.nchannels, PURPLE);
		}
	}
	for (const auto &w : map->walls) {
		draw_line(w.P0.x, w.P0.y, w.P1.x, w.P1.y, image.data, image.width, image.height, image.nchannels, PURPLE);

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

	Sitemap sitemap = {dist(gen), SITE_AREA};

	print_site(&sitemap);

	return 0;
}
