#include <iostream>
#include <random>
#include <queue>
#include <algorithm>
#include <unordered_map>
#include <list>
#include <functional>
#include <math.h>
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

glm::vec2 polygon_centroid(std::vector<glm::vec2> &vertices)
{
	glm::vec2 centroid = {0.f, 0.f};
	float area = 0.f;
	float x0 = 0.f; // Current vertex X
	float y0 = 0.f; // Current vertex Y
	float x1 = 0.f; // Next vertex X
	float y1 = 0.f; // Next vertex Y
	float a = 0.f;  // Partial signed area

	// For all vertices except last
	int i = 0;
	for (i = 0; i < vertices.size()-1; ++i) {
		x0 = vertices[i].x;
		y0 = vertices[i].y;
		x1 = vertices[i+1].x;
		y1 = vertices[i+1].y;
		a = x0*y1 - x1*y0;
		area += a;
		centroid.x += (x0 + x1)*a;
		centroid.y += (y0 + y1)*a;
	}

	// Do last vertex separately to avoid performing an expensive
	// modulus operation in each iteration.
	x0 = vertices[i].x;
	y0 = vertices[i].y;
	x1 = vertices[0].x;
	y1 = vertices[0].y;
	a = x0*y1 - x1*y0;
	area += a;
	centroid.x += (x0 + x1)*a;
	centroid.y += (y0 + y1)*a;

	area *= 0.5f;
	centroid.x /= (6.f*area);
	centroid.y /= (6.f*area);

	return centroid;
}

static bool comp_area(const struct section *a, const struct section*b)
{
	return a->area > b->area;
}

float polygon_area(std::vector<glm::vec2> &vertices)
{
	// triangle
	if (vertices.size() == 3) {
		return triangle_area(vertices[0], vertices[1], vertices[2]);
	}

	float area = 0.f;

	// Calculate value of shoelace formula
	int j = vertices.size() - 1;
	for (int i = 0; i < vertices.size(); i++) {
		area += (vertices[j].x + vertices[i].x) * (vertices[j].y - vertices[i].y);
		j = i;
	}

	return fabs(area / 2.f);
}

static bool winding(struct junction *a, struct junction *b, glm::vec2 center)
{
	glm::vec2 av = glm::normalize(a->position - center);
	glm::vec2 ab = glm::normalize(b->position - center);

	return atan2(av.y, av.x) < atan2(ab.y, ab.x);
}

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
			.street = false
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
		// sort the polygon vertices in order
		// ugly syntax
		std::sort(d.junctions.begin(), d.junctions.end(), std::bind(winding, std::placeholders::_1, std::placeholders::_2, d.center));
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
			std::queue<struct junction*> queue;
			queue.push(outward);
			glm::vec2 dir = glm::normalize(outward->position - inward->position);
			while (!queue.empty()) {
				struct junction *node = queue.front();
				queue.pop();
				struct junction *next = nullptr;
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
					node->street = true;
					next->street = true;
				}
			}
		}
	}
}

void Sitemap::divide_parcels(void)
{
	for (auto &d : districts) {
		if (d.radius > 0 && d.radius < TOWN_DISTRICT_RADIUS) {
			std::list<glm::vec2> polygon;
			for (auto jun : d.junctions) {
				printf("%f, %f\n", jun->position.x, jun->position.y);
				glm::vec2 p = 0.9f * (jun->position - d.center);
				polygon.push_back(d.center+p);
			}
			// set density
			for (std::list<glm::vec2>::iterator it = polygon.begin(); it != polygon.end();) {
				std::list<glm::vec2>::iterator next = std::next(it);
				if (next == polygon.end()) {
					next = polygon.begin();
				}
				glm::vec2 a = *it;
				glm::vec2 b = *next;
				//printf("distance %f\n", glm::distance(a, b));
				if (glm::distance(a, b) > 25.f) {
					polygon.insert(next, segment_midpoint(a, b));
				}
				if (next != polygon.begin()) {
					it = next;
				} else {
					it = polygon.end();
				}
			}
			divide_polygons(polygon, &d);
		}
	}
}

struct chainsplit {
	std::list<glm::vec2>::iterator target;
	std::list<glm::vec2>::iterator a;
	std::list<glm::vec2>::iterator b;
	glm::vec2 point;
};
			
// ugly
static struct chainsplit find_chainsplit(std::list<glm::vec2> &polygon)
{
	struct chainsplit split;

	// sort from longest to shortest segment
	float maxdistance = std::numeric_limits<float>::min();
	for (std::list<glm::vec2>::iterator it = polygon.begin(); it != polygon.end(); ++it) {
		std::list<glm::vec2>::iterator next = std::next(it);
		if (next == polygon.end()) {
			next = polygon.begin();
		}
		glm::vec2 a = *it;
		glm::vec2 b = *next;
		float dist = glm::distance(a, b);
		if (dist > maxdistance) {
			maxdistance = dist;
			split.a = it;
			split.b = next;
		}
	}

	// find the best possible way to split the polygon in two roughly equal areas
	float min = std::numeric_limits<float>::max();
	for (std::list<glm::vec2>::iterator it = polygon.begin(); it != polygon.end(); ++it) {
		if (it != split.a && it != split.b) {
			glm::vec2 point = *it;
			// now add the new vertex to the list
			const glm::vec2 splitpoint = closest_point_segment(point, *split.a, *split.b);
			const std::list<glm::vec2>::iterator splitstart = polygon.insert(split.b, splitpoint);
			// now split the polygon in two
			// first half
			std::list<glm::vec2> right;
			std::list<glm::vec2>::iterator chain = splitstart;
			while (chain != it) {
				right.push_back(*chain);
				chain++;	
				if (chain == polygon.end()) { 
					chain = polygon.begin(); 
				}
			}
			right.push_back(*it);
			// second half
			std::list<glm::vec2> left;
			chain = splitstart;
			while (chain != it) {
				left.push_back(*chain);
				if (chain == polygon.begin()) { 
					chain = std::prev(polygon.end()); 
				} else {
					chain--;	
				}
			}
			left.push_back(*it);

			//printf("right size %d\n", right.size());
			//printf("left size %d\n", left.size());

			std::vector<glm::vec2> rightpoints;
			std::vector<glm::vec2> leftpoints;
			for (auto point : right) {
				rightpoints.push_back(point);
			}
			for (auto point : left) {
				leftpoints.push_back(point);
			}
			float arearight = polygon_area(rightpoints);
			float arealeft = polygon_area(leftpoints);
			float ratio = arearight / arealeft;

			printf("area right %f\n", arearight);
			printf("area left %f\n", arealeft);
			printf("area ratio %f\n", ratio);

			// best ratio is the closest to 1
			float dist = fabs(1.F - ratio);
			if (dist < min) {
				min = dist;
				split.target = it;
				split.point = splitpoint;
			}

			// remove the inserted point from the original polygon
			polygon.erase(splitstart);
		}
	}

	return split;
}

static struct parcel make_parcel(glm::vec2 a, glm::vec2 b, glm::vec2 c, glm::vec2 d, const struct district *cell)
{
	struct parcel par;

	// find the front
	//for (const auto sect : cell->sections) {
	//}
	// find the back
	//
	par.frontleft = d;
	par.frontright = c;
	par.backleft = b;
	par.backright = a;
	par.owner = cell;

	std::vector<glm::vec2> vertices;
	vertices.push_back(par.frontleft);
	vertices.push_back(par.frontright);
	vertices.push_back(par.backleft);
	vertices.push_back(par.backright);
	par.centroid = polygon_centroid(vertices);
	par.direction = glm::normalize(segment_midpoint(par.frontleft, par.frontright) - par.centroid);
	return par;
}

void Sitemap::divide_polygons(std::list<glm::vec2> start, const struct district *cell)
{
	std::queue<std::list<glm::vec2>> queue; // queue of polygons to split
	queue.push(start);

	while (!queue.empty()) {
		std::list<glm::vec2> polygon = queue.front();
		queue.pop();
		if (polygon.size() == 4) {
			glm::vec2 a = polygon.front();
			polygon.pop_front();
			glm::vec2 b = polygon.front();
			polygon.pop_front();
			glm::vec2 c= polygon.front();
			polygon.pop_front();
			glm::vec2 d = polygon.front();
			polygon.pop_front();
			struct parcel par = make_parcel(a, b, c, d, cell);
			parcels.push_back(par);
		} else if (polygon.size() == 3) {
			//struct parcel par;
			//parcels.push_back(par);
		} else {
			// of all the polygon segments find a point that lies on a perpendicular line to the segment that divides the polygon in two "almost equal" areas 
			printf("\n");
			printf("%d polygon before\n", polygon.size());
			struct chainsplit split = find_chainsplit(polygon);
			printf("%d polygon after\n", polygon.size());
			// divide the polygon in two
			const std::list<glm::vec2>::iterator splitstart = polygon.insert(split.b, split.point);
			// first half
			std::list<glm::vec2> right;
			std::list<glm::vec2>::iterator it = splitstart;
			while (it != split.target) {
				right.push_back(*it);
				it++;	
				if (it == polygon.end()) { 
					it = polygon.begin(); 
				}
			}
			right.push_back(*split.target);
			// second half
			std::list<glm::vec2> left;
			it = splitstart;
			while (it != split.target) {
				left.push_back(*it);
				if (it == polygon.begin()) { 
					it = std::prev(polygon.end()); 
				} else {
					it--;	
				}
			}
			left.push_back(*split.target);

			std::vector<glm::vec2> rightpoints;
			std::vector<glm::vec2> leftpoints;
			for (auto point : right) {
				rightpoints.push_back(point);
			}
			for (auto point : left) {
				leftpoints.push_back(point);
			}
			float arearight = polygon_area(rightpoints);
			float arealeft = polygon_area(leftpoints);

			printf("area right %f\n", arearight);
			printf("area left %f\n", arealeft);

			if (arearight > 0.001f && arealeft > 0.001f) {
				queue.push(right);
				queue.push(left);
			}
			printf("right size %d\n", right.size());
			printf("left size %d\n", left.size());
		}
	}
}
