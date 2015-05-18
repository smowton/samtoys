targets: seektest subset bucket bamcmp rename_chroms reorder_chroms remove_qname_suffix filter_match_ratio

%: %.cpp
	g++ $^ -O3 -o $@ -std=c++11 -ggdb3 -lhts -lpthread
