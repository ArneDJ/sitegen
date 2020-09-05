CXX = g++

OBJDIR = obj
VPATH = src:

OBJECTS = main.o geom.o imp.o voronoi.o sitemap.o
OBJ = $(patsubst %,$(OBJDIR)/%,$(OBJECTS))

$(OBJDIR)/%.o: %.cpp
	$(CXX) -c -o $@ $<

site.out: $(OBJ)
	$(CXX) -o $@ $^

.PHONY : clean
clean:
	rm site.out -f $(OBJDIR)/*.o
