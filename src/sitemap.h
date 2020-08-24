enum SITEMAP_TYPE { VILLAGE, CASTLE, TOWN };

struct district;
struct junction;
struct section;

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
};

class Sitemap {
public:
	std::vector<struct district> districts;
	std::vector<struct junction> junctions;
	std::vector<struct section> sections;
	struct district *core;
	std::vector<struct segment> highways;
public:
	Sitemap(long seed, struct rectangle area);
	void make_diagram(void);
	void make_districts(void);
	void find_junction_radius(void);
	void outline_walls(void);
	void make_gateways(void);
	void make_highways(void);
private:
	long seed;
	struct rectangle area;
};
