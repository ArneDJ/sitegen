#include <iostream>
#include <random>
#include <queue>
#include <unordered_map>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "extern/stb_image_write.h"

#include "geom.h"
#include "imp.h"
#include "voronoi.h"

const struct rectangle SITE_AREA = { .min = {0.F, 0.F}, .max = {1024.F, 1024.F} };
const glm::vec2 CENTER = {512.F, 512.F};

void print_diagram(const Voronoi *voronoi) 
{
	std::unordered_map<int, bool> visited;
	std::unordered_map<int, int> depth;
	int root = 0;
	float min = 1024.f;
	for (const auto &cell : voronoi->cells) {
		visited[cell.index] = false;
		depth[cell.index] = 0;
		float dist = glm::distance(CENTER, cell.center);
		if (dist < min) { 
			min = dist;
			root = cell.index;
		}
	}

	int max = 0;
	std::queue<const struct cell*> queue;
	queue.push(&voronoi->cells[root]);
	while (!queue.empty()) {
		const struct cell *node = queue.front();
		queue.pop();
		int layer = depth[node->index] + 1;
		if (layer > max) { max = layer; }
		for (auto neighbor : node->neighbors) {
			if (visited[neighbor->index] == false) {
				visited[neighbor->index] = true;
				depth[neighbor->index] = layer;
				queue.push(neighbor);
			}
		}
	}
	struct byteimage image = blank_byteimage(3, 1024, 1024);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dist(0.f, 1.f);

	unsigned char ora[] = {255, 0, 0};
	unsigned char color[3];
	for (const auto &cell : voronoi->cells) {
		glm::vec3 rgb = {dist(gen), dist(gen), dist(gen)};
		float dist = 1.f - (depth[cell.index] / float(max));
		color[0] = 255 * dist;
		color[1] = 255 * dist;
		color[2] = 255 * dist;
		if (cell.index == root) {
			color[0] = 255;
			color[1] = 255;
			color[2] = 255;
		}

		glm::vec2 a = {round(cell.center.x), round(cell.center.y)};
		for (const auto &e : cell.edges) {
			// round points to rasterize properly
			glm::vec2 b = {round(e->v0->position.x), round(e->v0->position.y)};
			glm::vec2 c = {round(e->v1->position.x), round(e->v1->position.y)};
			draw_triangle(a, b, c, image.data, image.width, image.height, image.nchannels, color);
  		}

	}
	for (const auto &edg : voronoi->edges) {
		draw_line(edg.v0->position.x, edg.v0->position.y, edg.v1->position.x, edg.v1->position.y, image.data, image.width, image.height, image.nchannels, ora);
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
	std::uniform_real_distribution<float> dist_x(SITE_AREA.min.x, SITE_AREA.max.x);
	std::uniform_real_distribution<float> dist_y(SITE_AREA.min.y, SITE_AREA.max.y);

	std::vector<vec2> points;
	for (int i = 0; i < 100; i++) {
		points.push_back((vec2) {dist_x(gen), dist_y(gen)});
	}

	Voronoi voronoi;
	voronoi.gen_diagram(points, SITE_AREA.min, SITE_AREA.max, false);

	print_diagram(&voronoi);

	return 0;
}
