all: trilist

trilist: .obj/main.o .obj/Graph.o .obj/HashSet.o .obj/Timer.o .obj/Utility.o
	g++ .obj/main.o .obj/Graph.o .obj/HashSet.o .obj/Timer.o .obj/Utility.o -o trilist

.obj/main.o: main.cpp
	g++ -c -O3 -o .obj/main.o main.cpp

.obj/Graph.o: Graph.cpp
	g++ -I /usr/local/include/ -c -O3 -o .obj/Graph.o Graph.cpp

.obj/HashSet.o: HashSet.cpp
	g++ -c -O3 -o .obj/HashSet.o HashSet.cpp

.obj/Timer.o: Timer.cpp
	g++ -c -O3 -o .obj/Timer.o Timer.cpp
	
.obj/Utility.o: Utility.cpp
	g++ -c -O3 -o .obj/Utility.o Utility.cpp

clean:
	rm -rf *o .obj/
	mkdir .obj
