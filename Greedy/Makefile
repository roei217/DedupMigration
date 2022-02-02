CC = g++

# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS  = -std=c++11 -O3 -Wall
LINK = -lstdc++fs

GreedyLoadBalancerUnited: GreedyLoadBalancerUnited.o
	$(CC) $(CFLAGS) -o GreedyLoadBalancerUnited GreedyLoadBalancerUnited.o $(LINK)

clean:
	$(RM) GreedyLoadBalancerUnited.o GreedyLoadBalancerUnited
