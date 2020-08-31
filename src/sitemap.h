enum SITEMAP_TYPE { VILLAGE, CASTLE, TOWN };

struct district;
struct junction;
struct section;

struct parcel {
	// front faces the nearest street
	glm::vec2 frontleft;
	glm::vec2 frontright;
	// back
	glm::vec2 backleft;
	glm::vec2 backright;
	glm::vec2 centroid;
	glm::vec2 direction; // normalized direction vector to the closest street, if a parcel has a building on it it will be rotated in that direction
	const struct district *quarter;
};

struct section {
	int index;
	struct junction *j0 = nullptr;
	struct junction *j1 = nullptr;
	struct district *d0 = nullptr;
	struct district *d1 = nullptr;
	bool border;
	bool wall;
	float area;
	bool gateway;
};

struct junction {
	int index;
	glm::vec2 position;
	std::vector<struct junction*> adjacent;
	std::vector<struct district*> districts;
	bool border;
	int radius;
	bool wallcandidate;
	bool street;
};

struct district {
	int index;
	glm::vec2 center;
	std::vector<struct district*> neighbors;
	std::vector<struct junction*> junctions;
	std::vector<struct section*> sections;
	bool border;
	int radius; // distance to center in graph structure
	float area;
	bool wall;
	glm::vec2 centroid;
};

class Sitemap {
public:
	struct district *core;
	std::vector<struct segment> highways;
	std::vector<struct parcel> parcels;
	long seed;
	// graph structures
	std::vector<struct district> districts;
	std::vector<struct junction> junctions;
	std::vector<struct section> sections;
public:
	Sitemap(long seed, struct rectangle area);
	void make_diagram(void);
private:
	struct rectangle area;
private:
	void make_districts(void);
	void find_junction_radius(void);
	void outline_walls(void);
	void make_gateways(void);
	void make_highways(void);
	void divide_parcels(void);
	void divide_polygons(std::list<glm::vec2> start, const struct district *cell);
};
