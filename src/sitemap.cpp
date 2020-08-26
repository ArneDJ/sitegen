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
static const size_t MIN_GATEWAY_DISTANCE = 5;

Sitemap::Sitemap(long seed, struct rectangle area)
{
	this->seed = seed;
	this->area = area;

	make_diagram();

	make_districts();

	find_junction_radius();

	outline_walls();

	make_gateways();

	make_highways();

	divide_parcels();
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
		sections[index].gateway = false;
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

	/* TODO calculated area should be the rectangle INSIDE the walls */
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

static bool comp_area(const struct section *a, const struct section*b)
{
	return a->area > b->area;
}

void Sitemap::make_gateways(void)
{
	std::vector<struct section*> candidates;
	for (auto &sect : sections) {
		if (sect.wall) {
			candidates.push_back(&sect);
		}
	}
	// traverse per district
	std::unordered_map<int, bool> reserved;
	std::unordered_map<int, int> depth;
	for (const auto &dist : districts) {
		reserved[dist.index] = false;
		depth[dist.index] = 0;
	}
	// sort by largest area
	std::sort(candidates.begin(), candidates.end(), comp_area);
	for (auto sect : candidates) {
		if (reserved[sect->d0->index] == false && reserved[sect->d1->index] == false) {
			sect->gateway = true;
			std::queue<struct district*> queue;
			queue.push(sect->d0);
			queue.push(sect->d1);
			while (!queue.empty()) {
				struct district *node = queue.front();
				queue.pop();
				int layer = depth[node->index] + 1;
				reserved[node->index] = true;
				for (auto neighbor : node->neighbors) {
					if (reserved[neighbor->index] == false && layer < MIN_GATEWAY_DISTANCE) {
						depth[neighbor->index] = layer;
						queue.push(neighbor);
					}
				}
			}
		}
	}
}

void Sitemap::make_highways(void)
{
	// highway from gateway to town core
	for (auto &sect : sections) {
		if (sect.gateway) {
			highways.push_back((struct segment){sect.j0->position, sect.j1->position});
			struct junction *start = sect.j0->radius < sect.j1->radius ? sect.j0 : sect.j1;
			std::queue<struct junction*> queue;
			queue.push(start);
			while (!queue.empty()) {
				struct junction *node = queue.front();
				queue.pop();
				for (auto neighbor : node->adjacent) {
					if (neighbor->radius < node->radius) {
						queue.push(neighbor);
						highways.push_back((struct segment){node->position, neighbor->position});
						break;
					}
				}
			}
		}
	}

	// outside town center
	std::unordered_map<int, bool> visited;
	std::unordered_map<int, int> depth;
	for (const auto &j : junctions) {
		visited[j.index] = j.border;
		depth[j.index] = 0;
	}

	for (const auto &root : junctions) {
		if (root.border == true) {
			std::queue<const struct junction*> queue;
			queue.push(&root);
			while (!queue.empty()) {
				const struct junction *node = queue.front();
				queue.pop();
				int layer = depth[node->index] + 1;
				for (auto neighbor : node->adjacent) {
					if (visited[neighbor->index] == false) {
						visited[neighbor->index] = true;
						depth[neighbor->index] = layer;
						queue.push(neighbor);
					} else if (depth[neighbor->index] > layer) {
						depth[neighbor->index] = layer;
						queue.push(neighbor);
					}
				}
			}
		}
	}
	for (auto &sect : sections) {
		if (sect.gateway) {
			struct junction *outward = sect.j0->radius > sect.j1->radius ? sect.j0 : sect.j1;
			struct junction *inward = sect.j0->radius < sect.j1->radius ? sect.j0 : sect.j1;
			std::queue<const struct junction*> queue;
			queue.push(outward);
			glm::vec2 dir = glm::normalize(outward->position - inward->position);
			while (!queue.empty()) {
				const struct junction *node = queue.front();
				queue.pop();
				const struct junction *next = nullptr;
				float maxdot = -1.f;
				for (auto neighbor : node->adjacent) {
					if (depth[neighbor->index] < depth[node->index] && neighbor->wallcandidate == false) {
						glm::vec2 nextdir = glm::normalize(neighbor->position - node->position);
						float dotp = glm::dot(nextdir, dir);
						if (dotp > maxdot) {
							maxdot = dotp;
							next = neighbor;
						}
					}
				}
				if (next != nullptr) {
					queue.push(next);
					highways.push_back((struct segment) {node->position, next->position});
				}
			}
		}
	}
}

void Sitemap::divide_parcels(void)
{
	for (auto &d : districts) {
		if (d.radius < WALL_RADIUS) {
			std::list<glm::vec2> polygon;
			for (auto jun : d.junctions) {
				polygon.push_back(jun->position);
			}
			divide_polygons(polygon);
		}
	}
}

// a polygon can be stored as a double linked list of connected vertices
void Sitemap::divide_polygons(std::list<glm::vec2> start)
{
	std::queue<std::list<glm::vec2>> polygons; // queue of polygons to split
	polygons.push(start);

	while (!polygons.empty()) {
		std::list<glm::vec2> polygon = polygons.front();
		polygons.pop();
		if (polygon.size() > 4) {
			struct segment longest = {{0.f, 0.f} , {0.f, 0.f}};
			std::list<glm::vec2>::iterator longestv;
			float max = 0.F;
			// find the longest segment in the polygon
			for (std::list<glm::vec2>::iterator it = polygon.begin(); it != polygon.end(); ++it) {
				std::list<glm::vec2>::iterator next = std::next(it);
				if (next == polygon.end()) {
					next = polygon.begin();
				}
				glm::vec2 a = *it;
				glm::vec2 b = *next;
				float dist = glm::distance(a, b);
				if (dist > max) {
					max = dist;
					longest.P0 = a;
					longest.P1 = b;
					longestv = next;
				}
			}
			// find the point with the longest distance to the segment
			std::list<glm::vec2>::iterator target;
			max = 0.F;
			for (std::list<glm::vec2>::iterator it = polygon.begin(); it != polygon.end(); ++it) {
				glm::vec2 point = *it;
				float dist = sqrdist_point_segment(point, longest.P0, longest.P1);
				if (dist > max) {
					max = dist;
					target = it;
				}
			}
			// now add the new vertex to the list
			const glm::vec2 splitpoint = closest_point_segment(*target, longest.P0, longest.P1);
			const std::list<glm::vec2>::iterator splitstart = polygon.insert(longestv, splitpoint);
			// now split the polygon in two
			// first half
			std::list<glm::vec2> ying;
			std::list<glm::vec2>::iterator it = splitstart;
			while (it != target) {
				ying.push_back(*it);
				if (it == polygon.end()) { 
					it = polygon.begin(); 
				} else {
					it++;	
				}
			}
			// second half
			std::list<glm::vec2> yang;
			it = splitstart;
			while (it != target) {
				yang.push_back(*it);
				if (it == polygon.begin()) { 
					it = polygon.end(); 
				} else {
					it--;	
				}
			}
		printf("original polygon\n");
		for (auto &p : polygon) {
			printf("%f, %f\n", p.x, p.y);
		}
		printf("first half\n");
		for (auto &p : ying) {
			printf("%f, %f\n", p.x, p.y);
		}
		printf("second half\n");
		for (auto &p : yang) {
			printf("%f, %f\n", p.x, p.y);
		}
		} else {
			struct parcel par;
			par.a = polygon.front();
			polygon.pop_front();
			par.b = polygon.front();
			polygon.pop_front();
			par.c = polygon.front();
			polygon.pop_front();
			par.d = polygon.front();
			polygon.pop_front();
			parcels.push_back(par);
		}
	}
}
