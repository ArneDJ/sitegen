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

static const size_t MIN_TOWNGATE_DISTANCE = 6; // minimum distance of wall segments between two town gates 
static const float MIN_GATEWALL_DISTANCE = 90.F;
static const float MIN_INWARD_AREA = 1500.F;

// checks how "isoscelene" a triangle is by comparing the difference in length of two triangle legs, the smaller the number the more "isoscelene" it is
static float triangle_isoscelenes(glm::vec2 a, glm::vec2 b, glm::vec2 c)
{
	float d0 = glm::distance(a, b);
	float d1 = glm::distance(b, c);

	float circumfence = d0 + glm::distance(b, c) + d1;
	return fabs(d0 - d1) / circumfence;
}

float triangle_area(glm::vec2 a, glm::vec2 b, glm::vec2 c)
{
	float area = ((b.x - a.x)*(c.y - a.y) - (c.x - a.x)*(b.y - a.y))/2.0;
	return (area > 0.0) ? area : -area;
}

struct meta {
	const struct junction *j;
	const struct junction *prev;
	const struct junction *next;
	float iso;
};

bool lowest(struct meta &a, struct meta &b) 
{
	return a.iso < b.iso;
}

bool highest_area(struct meta &a, struct meta &b) 
{
	return triangle_area(a.j->position, a.prev->position, a.next->position) > triangle_area(b.j->position, b.prev->position, b.next->position);
}

static struct towngate make_gate(const struct meta *data)
{
	struct towngate gate;
	gate.inward = data->j;
	gate.outward = nullptr;

	int max = 0;
	struct district *cell;
	for (const auto d : data->j->districts) {
		if (d->radius > max) {
			max = d->radius;
			cell = d;
		}
	}

	gate.wall = {data->prev->position, data->next->position};
	gate.gatepoint = segment_midpoint(gate.wall.P0, gate.wall.P1);

	float mindot = 1.f;
	glm::vec2 a = glm::normalize(gate.wall.P0 - gate.gatepoint);
	for (const auto j : cell->junctions) {
		if (j != data->j && j != data->prev && j != data->next && j->wallcandidate == false) {
			glm::vec2 b = glm::normalize(j->position - gate.gatepoint);
			float dotp = glm::dot(a, b);
			float dist = fabs(dotp);
			if (dist < mindot && dist < 0.5F) {
				mindot = dist;
				gate.outward = j;
			}
		}
	}

	return gate;
}

Sitemap::Sitemap(long seed, struct rectangle area)
{
	this->seed = seed;
	this->area = area;

	adapt_diagram();

	std::unordered_map<int, bool> visited;
	//struct district *centerdistrict;
	float min = glm::distance(this->area.min, this->area.max);
	glm::vec2 center = 0.5f * this->area.max;

	for (auto &d : districts) {
		d.radius = 0;
		visited[d.index] = false;
		float dist = glm::distance(center, d.center);
		if (dist < min) {
			min = dist;
			towncenter = &d;
		}
	}

	std::queue<const struct district*> queue;
	queue.push(towncenter);
	visited[towncenter->index] = true;
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

	for (auto &j : junctions) {
		j.radius = 0;
		for (const auto &d : j.districts) {
			j.radius += d->radius;
		}
		if (j.radius > 6 && j.radius < 9 && j.border == false) {
			j.wallcandidate = true;
		} else {
			j.wallcandidate = false;
		}
	}

	outline_walls();

	make_highways();
}

void Sitemap::adapt_diagram(void)
{
	std::mt19937 gen(seed);
	std::uniform_real_distribution<float> dist_x(area.min.x, area.max.x);
	std::uniform_real_distribution<float> dist_y(area.min.y, area.max.y);

	std::vector<glm::vec2> points;
	for (int i = 0; i < 128; i++) {
		points.push_back((glm::vec2) {dist_x(gen), dist_y(gen)});
	}

	Voronoi voronoi;
	voronoi.gen_diagram(points, area.min, area.max, true);

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
			.border = false
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
}

void Sitemap::outline_walls(void)
{
	// make the walls
	std::unordered_map<int, bool> discovered;
	std::list<struct junction*> chain;
	for (const auto &j : junctions) {
		discovered[j.index] = false;
	}
	for (auto &sect : sections) {
		if (sect.j0->wallcandidate == true && sect.j1->wallcandidate == true && discovered[sect.j0->index] == false && discovered[sect.j1->index] == false) {
			discovered[sect.j0->index] = true;
			discovered[sect.j1->index] = true;
			std::list<struct junction*> queue;
			chain.push_back(sect.j0);
			chain.push_front(sect.j1);
			queue.push_back(sect.j0);
			while (!queue.empty()) {
				struct junction *node = queue.front();
				queue.pop_front();
				for (auto neighbor : node->adjacent) {
					if (discovered[neighbor->index] == false && neighbor->wallcandidate == true) {
						discovered[neighbor->index] = true;
						queue.push_back(neighbor);
						chain.push_back(neighbor);
					}
				}
			}
		}
	}

	std::vector<struct meta> data;
	for (std::list<struct junction*>::iterator it = chain.begin(); it != chain.end(); ++it) {
		auto nx = std::next(it);
		if (nx == chain.end()) {
			nx = chain.begin();
		}
		auto nxnx = std::next(nx);
		if (nxnx == chain.end()) {
			nxnx = chain.begin();
		}
		const struct junction *a = *it;
		const struct junction *b = *nx;
		const struct junction *c = *nxnx;
		float iso = triangle_isoscelenes(a->position, b->position, c->position);
		// check if the triangle points "inside"
		if (b->radius < 8) {
			//if (iso < 0.5F) {
				data.push_back((struct meta) {b, a, c, iso});
			//}
			it = chain.erase(nx);
			it = std::prev(it);
		}
	}

	// make town gates
	// reset discovered
	std::unordered_map<int, int> depth;
	for (const auto &j : junctions) {
		discovered[j.index] = false;
		depth[j.index] = 0;
	}
	std::sort(data.begin(), data.end(), highest_area);
	for (auto &d : data) {
		if (discovered[d.j->index] == false) {
			discovered[d.j->index] = true;
			entrances.push_back(d.j);
			struct towngate gate = make_gate(&d);
			float area = triangle_area(gate.wall.P0, gate.wall.P1, gate.inward->position);
			if (area > MIN_INWARD_AREA && gate.outward != nullptr) {
				towngates.push_back(gate);
				std::queue<const struct junction*> queue;
				queue.push(d.j);
				while (!queue.empty()) {
					const struct junction *node = queue.front();
					queue.pop();
					int layer = depth[node->index] + 1;
					if (layer > MIN_TOWNGATE_DISTANCE) {
						break;
					}
					for (auto neighbor : node->adjacent) {
						if (neighbor->wallcandidate == true && discovered[neighbor->index] == false) {
							discovered[neighbor->index] = true;
							depth[neighbor->index] = layer;
							queue.push(neighbor);
						}
					}
				}
			}
		}
	}

	for (std::list<struct junction*>::iterator it = chain.begin(); it != chain.end(); ++it) {
		auto nx = std::next(it);
		if (nx == chain.end()) {
			nx = chain.begin();
		}

		const struct junction *j = *it;
		const struct junction *nxj = *nx;

		struct segment S = { j->position, nxj->position };
		walls.push_back(S);
	}
}

void Sitemap::make_highways(void)
{
	std::unordered_map<int, bool> visited;
	std::unordered_map<int, int> depth;
	for (const auto &j : junctions) {
		visited[j.index] = false;
		depth[j.index] = 0;
	}
	for (const auto &root : towncenter->junctions) {
		visited[root->index] = true;
	}

	// town center
	for (const auto &root : towncenter->junctions) {
		std::queue<const struct junction*> queue;
		queue.push(root);
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
	for (auto &gate : towngates) {
		std::queue<const struct junction*> queue;
		queue.push(gate.inward);
		glm::vec2 dir = glm::normalize(gate.inward->position - gate.gatepoint);
		while (!queue.empty()) {
			const struct junction *node = queue.front();
			queue.pop();
			const struct junction *next = nullptr;
			float maxdot = -1.f;
			for (auto neighbor : node->adjacent) {
				if (depth[neighbor->index] < depth[node->index]) {
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

	// outside center
	// reset
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
	for (auto &gate : towngates) {
		std::queue<const struct junction*> queue;
		queue.push(gate.outward);
		glm::vec2 dir = glm::normalize(gate.outward->position - gate.gatepoint);
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
