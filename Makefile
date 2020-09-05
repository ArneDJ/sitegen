CXX = g++
VPATH = src:
OBJECTS = main.o geom.o imp.o voronoi.o sitemap.o

site.out : $(OBJECTS)
	$(CXX) -o site.out $(OBJECTS)

main.o : main.cpp
geom.o : geom.cpp
imp.o : imp.cpp
voronoi.o : voronoi.cpp
sitemap.o : sitemap.cpp

.PHONY : clean
clean:
	rm site.out -f *.o
