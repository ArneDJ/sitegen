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

const struct rectangle SITE_AREA = { .min = {0.F, 0.F}, .max = {1024.F, 1024.F} };
const glm::vec2 CENTER = {512.F, 512.F};

bool city_edge(const struct vertex *v, std::unordered_map<int, int> &depth)
{
	bool inner = false;
	bool outer = false;
	for (auto &c : v->cells) {
		int d = depth[c->index];
		if (d == 2) {
			inner = true;
		}
		if (d == 3) {
			outer = true;
		}
	}

	return inner && outer;
}

void print_paths(const Voronoi *voronoi, std::vector<const struct vertex*> &gatepoints, std::unordered_map<int, int> &depth, const struct cell *center, struct byteimage *image)
{
	unsigned char pur[] = {255, 0, 255};
	unsigned char grn[] = {0, 255, 0};
	for (const auto start : gatepoints) {
		plot(start->position.x, start->position.y, image->data, image->width, image->height, image->nchannels, pur);
	}

	// construct town diagram excluding nodes on that are part of the wall
	std::vector<const struct vertex*> diagram;
	std::unordered_map<int, bool> visited;
	std::unordered_map<int, int> streetdepth;
	std::unordered_map<int, bool> forbidden;
	std::unordered_map<int, int> distcenter;
	for (auto &v : voronoi->vertices) {
		int total = 0;
		for (auto c : v.cells) {
			const int d = depth[c->index];
			total += d;
		}
		distcenter[v.index] = total;
		if (total < 8) {
			diagram.push_back(&v);
			visited[v.index] = false;
			streetdepth[v.index] = 0;
			forbidden[v.index] = false;
		} else {
			forbidden[v.index] = true;
		}
	}

	for (auto root : center->vertices) {
		forbidden[root->index] = true;
	}
	// now do breath first search from the town center
	for (auto root : center->vertices) {
		printf("%d\n", distcenter[root->index]);
		plot(root->position.x, root->position.y, image->data, image->width, image->height, image->nchannels, pur);
		std::queue<const struct vertex*> queue;
		queue.push(root);
		while (!queue.empty()) {
			const struct vertex *node = queue.front();
			queue.pop();
			int layer = streetdepth[node->index] + 1;
			for (auto neighbor : node->adjacent) {
				if (forbidden[neighbor->index] == false) {
					if (visited[neighbor->index] == false) {
						streetdepth[neighbor->index] = layer;
						visited[neighbor->index] = true;
						queue.push(neighbor);
					} else if (streetdepth[neighbor->index] > layer) {
						streetdepth[neighbor->index] = layer;
						queue.push(neighbor);
					}
				}
			}
		}
	}
	for (auto &v : voronoi->vertices) {
		visited[v.index] = false;
	}
	for (const auto start : gatepoints) {
		std::queue<const struct vertex*> queue;
		queue.push(start);
		while (!queue.empty()) {
			const struct vertex *node = queue.front();
			queue.pop();
			int layer = streetdepth[node->index];
			for (auto neighbor : node->adjacent) {
				if (forbidden[neighbor->index] == false) {
					//printf("%d\n", streetdepth[neighbor->index]);
					if (streetdepth[neighbor->index] <= layer && visited[neighbor->index] == false) {
						visited[neighbor->index] = true;
						draw_line(node->position.x, node->position.y, neighbor->position.x, neighbor->position.y, image->data, image->width, image->height, image->nchannels, grn);
						queue.push(neighbor);
						break;
					}
				} else if (distcenter[neighbor->index] < 3) {
						draw_line(node->position.x, node->position.y, neighbor->position.x, neighbor->position.y, image->data, image->width, image->height, image->nchannels, grn);
				}
			}
		}
	}
}

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
	visited[voronoi->cells[root].index] = true;
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
	printf("%d\n", depth[voronoi->cells[root].index]);

	std::unordered_map<int, bool> discovered;
	std::list<const struct vertex*> chain;
	for (auto &v : voronoi->vertices) {
		discovered[v.index] = false;
	}
	for (auto &edg : voronoi->edges) {
		if (edg.c0 != nullptr && edg.c1 != nullptr) {
			int d0 = depth[edg.c0->index];
			int d1 = depth[edg.c1->index];
			if ((d0 == 2 && d1 == 3) || (d1 == 2 && d0 == 3)) {
				if (discovered[edg.v0->index] == false && discovered[edg.v1->index] == false) {
				discovered[edg.v0->index] = true;
				discovered[edg.v1->index] = true;
				chain.push_back(edg.v0);
				//chain.push_back(edg.v1);
				std::list<const struct vertex*> queue;
				queue.push_back(edg.v1);
				while (!queue.empty()) {
					const struct vertex *node = queue.front();
					chain.push_back(node);
					queue.pop_front();
					for (auto neighbor : node->adjacent) {
						if (discovered[neighbor->index] == false && city_edge(neighbor, depth) == true) {
							discovered[neighbor->index] = true;
							queue.push_back(neighbor);
						}
					}
				}
				}
			}
		}
	}

	struct byteimage image = blank_byteimage(3, 1024, 1024);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dist(0.f, 1.f);

	unsigned char ora[] = {255, 0, 0};
	unsigned char blu[] = {0, 0, 255};
	unsigned char pur[] = {255, 0, 255};
	unsigned char grn[] = {0, 255, 0};
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

	int count = 0;
	std::vector<glm::vec2> gate_points;
	std::vector<const struct vertex*> gate_street;
	for (std::list<const struct vertex*>::iterator it = chain.begin(); it != chain.end(); ++it) {
		auto nx = std::next(it);
		if (nx == chain.end()) {
			nx = chain.begin();
		}
		auto nxnx = std::next(nx);
		const struct vertex *a = *it;
		const struct vertex *b = *nx;
		const struct vertex *c = *nxnx;
		// check if the triangel points "inside"
		int total = 0;
		for (auto c : b->cells) {
			const int d = depth[c->index];
			total += d;
		}
		if (total < 8) {
			it = chain.erase(nx);
			if (count > 5) {
				glm::vec2 mid = segment_midpoint(a->position, c->position);
				gate_points.push_back(mid);
				gate_street.push_back(b);
				count = 0;
			}
		}
		count++;
	}
	for (std::list<const struct vertex*>::iterator it = chain.begin(); it != chain.end(); ++it) {
		auto nx = std::next(it);
		if (nx == chain.end()) {
			nx = chain.begin();
		}

		const struct vertex *v = *it;
		const struct vertex *nxv = *nx;

		draw_line(v->position.x, v->position.y, nxv->position.x, nxv->position.y, image.data, image.width, image.height, image.nchannels, blu);
	}
	for (const auto &mid : gate_points) {
		plot(mid.x, mid.y, image.data, image.width, image.height, image.nchannels, grn);
		plot(mid.x+1, mid.y, image.data, image.width, image.height, image.nchannels, grn);
		plot(mid.x-1, mid.y, image.data, image.width, image.height, image.nchannels, grn);
		plot(mid.x, mid.y+1, image.data, image.width, image.height, image.nchannels, grn);
		plot(mid.x, mid.y-1, image.data, image.width, image.height, image.nchannels, grn);
	}
	
	print_paths(voronoi, gate_street, depth, &voronoi->cells[root], &image);

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
	voronoi.gen_diagram(points, SITE_AREA.min, SITE_AREA.max, true);

	print_diagram(&voronoi);

	return 0;
}
