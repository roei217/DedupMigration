CC = g++
OBJS = Main.o AlgorithmDSManager.o HierarchicalClustering.o Node.o Utility.o CommandLineParser.o
OBJSGCH = AlgorithmDSManager.hpp.gch HierarchicalClustering.hpp.gch Node.hpp.gch Utility.hpp.gch CommandLineParser.hpp.gch
EXEC = hc
DEBUG_FLAG = #do -g for debugging
LINKING = #
CFLAGS = -std=c++14 $(LINKING) -O3

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(DEBUG_FLAG) $(OBJS) -o $@
	
Main.o: Main.cpp
	$(CC) $(CFLAGS) -c Main.cpp
AlgorithmDSManager.o: AlgorithmDSManager.cpp AlgorithmDSManager.hpp
	$(CC) $(CFLAGS) -c AlgorithmDSManager.cpp
HierarchicalClustering.o: HierarchicalClustering.cpp HierarchicalClustering.hpp
	$(CC) $(CFLAGS) -c HierarchicalClustering.cpp
Node.o: Node.cpp Node.hpp
	$(CC) $(CFLAGS) -c Node.cpp
Utility.o: Utility.cpp Utility.hpp
	$(CC) $(CFLAGS) -c Utility.cpp
CommandLineParser.o: CommandLineParser.cpp CommandLineParser.hpp
	$(CC) $(CFLAGS) -c CommandLineParser.cpp

clean:
	rm -f $(OBJS)
	rm -f $(OBJSGCH)