#include <iostream>
#include <random>
#include <queue>
#include <algorithm>
#include <unordered_map>
#include <list>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>

#include "geom.h"
#include "imp.h"
#include "voronoi.h"
#include "sitemap.h"

static const size_t MAX_CELLS = 128;
static const size_t VILLAGE_DISTRICT_RADIUS = 2;
static const size_t TOWN_DISTRICT_RADIUS = 3;
static const size_t FIELD_DISTRICT_RADIUS = 5;
static const size_t WALL_RADIUS = 2;

Sitemap::Sitemap(long seed, struct rectangle area)
{
	this->seed = seed;
	this->area = area;

	make_diagram();

	make_districts();

	find_junction_radius();

	outline_walls();
}

void Sitemap::make_diagram(void)
{
	std::mt19937 gen(seed);
	std::uniform_real_distribution<float> dist_x(area.min.x, area.max.x);
	std::uniform_real_distribution<float> dist_y(area.min.y, area.max.y);

	std::vector<glm::vec2> locations;
	for (int i = 0; i < MAX_CELLS; i++) {
		locations.push_back((glm::vec2) {dist_x(gen), dist_y(gen)});
	}

	Voronoi voronoi;
	voronoi.gen_diagram(locations, area.min, area.max, true);

	districts.resize(voronoi.cells.size());
	junctions.resize(voronoi.vertices.size());
	sections.resize(voronoi.edges.size());

	// adopt cell structures
	for (const auto &cell : voronoi.cells) {
		std::vector<struct district*> dneighbors;
		for (const auto &neighbor : cell.neighbors) {
			dneighbors.push_back(&districts[neighbor->index]);
		}
		std::vector<struct junction*> djunctions;
		for (const auto &vertex : cell.vertices) {
			djunctions.push_back(&junctions[vertex->index]);
		}
		std::vector<struct section*> dsections;
		for (const auto &edge : cell.edges) {
			dsections.push_back(&sections[edge->index]);
		}

		struct district d = {
			.index = cell.index,
			.center = cell.center,
			.neighbors = dneighbors,
			.junctions = djunctions,
			.sections = dsections,
			.border = false,
			.radius = 0,
			.area = 0.f,
			.wall = false,
		};

		districts[cell.index] = d;
	}

	// adapt vertex structures
	for (const auto &vertex : voronoi.vertices) {
		std::vector<struct junction*> adjacent;
		for (const auto &neighbor : vertex.adjacent) {
			adjacent.push_back(&junctions[neighbor->index]);
		}
		std::vector<struct district*> touches;
		for (const auto &cell : vertex.cells) {
			touches.push_back(&districts[cell->index]);
		}

		struct junction c = {
			.index = vertex.index,
			.position = vertex.position,
			.adjacent = adjacent,
			.districts = touches,
			.border = false,
		};

		junctions[vertex.index] = c;
	}

	// adapt edge structures
	for (const auto &edge : voronoi.edges) {
		int index = edge.index;
		sections[index].j0 = &junctions[edge.v0->index];
		sections[index].j1 = &junctions[edge.v1->index];
		sections[index].border = false;
		sections[index].wall = false;
		sections[index].area = 0.f;
		if (edge.c0 != nullptr) {
			sections[index].d0 = &districts[edge.c0->index];
		} else {
			sections[index].d0 = &districts[edge.c1->index];
			sections[index].d0->border = true;
			sections[index].border = true;
			sections[index].j0->border = true;
			sections[index].j1->border = true;
		}
		if (edge.c1 != nullptr) {
			sections[index].d1 = &districts[edge.c1->index];
		} else {
			sections[index].d1 = &districts[edge.c0->index];
			sections[index].d1->border = true;
			sections[index].border = true;
			sections[index].j0->border = true;
			sections[index].j1->border = true;
		}
	}

	for (auto &d : districts) {
		d. area = 0.f;
		for (const auto &s : d.sections) {
			d.area += triangle_area(d.center, s->j0->position, s->j1->position);
		}
	}
}

void Sitemap::make_districts(void)
{
	// find the center tile
	std::unordered_map<int, bool> visited;
	float min = glm::distance(area.min, area.max);
	glm::vec2 center = 0.5f * area.max;

	for (auto &d : districts) {
		visited[d.index] = false;
		d.radius = 0;
		float dist = glm::distance(center, d.center);
		if (dist < min) {
			min = dist;
			core = &d;
		}
	}

	// distance to center tile in graph
	std::queue<const struct district*> queue;
	queue.push(core);
	visited[core->index] = true;
	while (!queue.empty()) {
		const struct district *node = queue.front();
		queue.pop();
		int depth = node->radius + 1;
		for (auto neighbor : node->neighbors) {
			if (visited[neighbor->index] == false) {
				visited[neighbor->index] = true;
				neighbor->radius = depth;
				queue.push(neighbor);
			}
		}
	}
}

void Sitemap::find_junction_radius(void)
{
	std::unordered_map<int, bool> visited;
	for (auto &j : junctions) {
		visited[j.index] = false;
		j.radius = 0;
	}
	for (const auto &root : core->junctions) {
		visited[root->index] = true;
	}

	// town center
	for (const auto &root : core->junctions) {
		std::queue<const struct junction*> queue;
		queue.push(root);
		while (!queue.empty()) {
			const struct junction *node = queue.front();
			queue.pop();
			int layer = node->radius + 1;
			for (auto neighbor : node->adjacent) {
				if (visited[neighbor->index] == false) {
					visited[neighbor->index] = true;
					neighbor->radius = layer;
					queue.push(neighbor);
				} else if (neighbor->radius > layer) {
					neighbor->radius = layer;
					queue.push(neighbor);
				}
			}
		}
	}
}

void Sitemap::outline_walls(void)
{
	for (auto &sect : sections) {
		if (sect.d0->radius == WALL_RADIUS && sect.d1->radius == WALL_RADIUS) {
			sect.wall = true;
			sect.d0->wall = true;
			sect.d1->wall = true;
		} else {
			sect.wall = false;
		}
	}

	for (auto &sect : sections) {
		sect.area = 0.f;
		if (sect.wall) {
			glm::vec2 outward = segment_midpoint(sect.d0->center, sect.d1->center);
			glm::vec2 right = sect.d0->center - outward;
			glm::vec2 a = sect.j0->position - right;;
			glm::vec2 b = sect.j0->position + right;;
			glm::vec2 c = sect.j1->position - right;;
			glm::vec2 d = sect.j1->position + right;;
			sect.area += triangle_area(a, c, d);
			sect.area += triangle_area(a, d, b);
		}
	}

	for (auto &j : junctions) {
		// look if all 3 of the districts have walls
		bool walls = true;
		for (auto &d : j.districts) {
			if (d->wall == false) {
				walls = false;
			}
		}
		// if it does remove the wall that is the closest to the core
		if (walls) {
			int minradius = MAX_CELLS;
			struct district *target = nullptr;
			for (auto d : j.districts) {
				int max = 0;
				for (auto &neighbor : d->neighbors) {
					if (neighbor->radius > max) {
						max = neighbor->radius;
					}
				}
				if (max < minradius) {
					minradius = max;
					target = d;
				}
			}
			if (target != nullptr) {
				target->wall = false;
				for (auto targetsect : target->sections) {
					targetsect->wall = false;
				}
			}
		}
	}
}
