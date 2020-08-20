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
};

struct junction {
	int index;
	glm::vec2 position;
	std::vector<struct junction*> adjacent;
	std::vector<struct district*> districts;
	bool border;
	int radius;
};

struct district {
	int index;
	glm::vec2 center;
	std::vector<struct district*> neighbors;
	std::vector<struct junction*> junctions;
	std::vector<struct section*> sections;
	bool border;
	int radius; // distance to center in graph structure
};

class Sitemap {
public:
	std::vector<struct district> districts;
	std::vector<struct junction> junctions;
	std::vector<struct section> sections;
	std::vector<struct segment> walls;
public:
	Sitemap(long seed, struct rectangle area);
	void adapt_diagram(void);
private:
	long seed;
	struct rectangle area;
};
