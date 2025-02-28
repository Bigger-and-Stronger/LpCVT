CXX=g++
COPT=-I../ -frounding-math -O3 -DCVT_MULTITHREAD 
OBJS=processor.o delaunay_CGAL.o delaunay.o mesh.o RVD_predicates.o F_Lp.o main.o
LIBS=-lCGAL 
LDXX=g++

LpCVT: $(OBJS)
	$(LDXX) $(LDOPTS) $(LIBS) $(OBJS) -o LpCVT

clean:
	rm -f $(OBJS) LpCVT *.obj *.eobj

processor.o: common/processor.cpp 
	$(CXX) -c $(COPT) common/processor.cpp

delaunay_CGAL.o: combinatorics/delaunay_CGAL.cpp 
	$(CXX) -c $(COPT) combinatorics/delaunay_CGAL.cpp 

delaunay.o: combinatorics/delaunay.cpp 
	$(CXX) -c $(COPT) combinatorics/delaunay.cpp 

mesh.o: combinatorics/mesh.cpp 
	$(CXX) -c $(COPT) combinatorics/mesh.cpp 

RVD_predicates.o: combinatorics/exact/RVD_predicates.cpp
	$(CXX) -c $(COPT) combinatorics/exact/RVD_predicates.cpp

F_Lp.o: algebra/F_Lp.cpp
	$(CXX) -c $(COPT) algebra/F_Lp.cpp

main.o: main.cpp
	$(CXX) -c $(COPT) main.cpp
