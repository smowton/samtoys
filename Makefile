targets: seektest subset bucket

%: %.cpp
	g++ $^ -O3 -o $@ -std=c++11 -ggdb3 -lhts -lpthread
