

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <vector>
#include <algorithm>
#include <map>

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

// Borrowed from Samtools source, since samtools sort -n uses this ordering:

static int strnum_cmp(const char *_a, const char *_b)
{
    const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
    const unsigned char *pa = a, *pb = b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            while (*pa == '0') ++pa;
            while (*pb == '0') ++pb;
            while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
            if (isdigit(*pa) && isdigit(*pb)) {
                int i = 0;
                while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
                return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
            } else if (isdigit(*pa)) return 1;
            else if (isdigit(*pb)) return -1;
            else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
        } else {
            if (*pa != *pb) return (int)*pa - (int)*pb;
            ++pa; ++pb;
        }
    }
    return *pa? 1 : *pb? -1 : 0;
}

static int flag2mate(const bam1_t* rec) {

  if(rec->core.flag & BAM_FREAD1)
    return 1;
  else if(rec->core.flag & BAM_FREAD2)
    return 2;
  return 0;

}

static bool mate_eq(const bam1_t* a, const bam1_t* b) {

  return flag2mate(a) == flag2mate(b);

}

static bool mate_lt(const bam1_t* a, const bam1_t* b) {

  return flag2mate(a) < flag2mate(b);

}

struct BamRecVector {

  std::vector<bam1_t*> recs;

  BamRecVector() {}
  ~BamRecVector() {
    clear();
  }

  void copy_add(bam1_t* src) {
    recs.push_back(bam_dup1(src));
  }

  void clear() {
    for(std::vector<bam1_t*>::iterator it = recs.begin(), itend = recs.end(); it != itend; ++it)
      bam_destroy1(*it);
    recs.clear();
  }

  void sort() {
    std::sort(recs.begin(), recs.end(), mate_lt);
  }

};

struct CmpIntVector {

  std::vector<int> data;

  bool operator==(const CmpIntVector& other) const {

    if(other.data.size() != data.size())
      return false;

    for(int i = 0, ilim = data.size(); i != ilim; ++i)
      if(data[i] != other.data[i])
	return false;

    return true;

  }

  bool operator<(const CmpIntVector& other) const {

    for(int i = 0, ilim = std::min(data.size(), other.data.size()); i != ilim; ++i) {

      if(data[i] < other.data[i])
	return true;
      else if(data[i] > other.data[i])
	return false;

    }

    return  other.data.size() > data.size();

  }

};

int main(int argc, char** argv) {
  
  if(argc < 3) {

    fprintf(stderr, "Usage: filter_hits in.xam out.txt\n");
    exit(1);

  }

  htsFile* hfi = hts_open(argv[1], "r");
  if(!hfi) {
    fprintf(stderr, "Failed to open %s\n", argv[1]);
    exit(1);
  }

  FILE* fo = fopen(argv[2], "w");
  if(!fo) {
    fprintf(stderr, "Failed to open %s\n", argv[2]);
    exit(1);
  }

  bam_hdr_t* header = sam_hdr_read(hfi);

  bam1_t* prevrec = bam_init1();
  bam1_t* rec = bam_init1();

  std::map<CmpIntVector, int> counts;

  CmpIntVector contigs;

  BamRecVector matchingRecs;

  do {

    if(sam_read1(hfi, header, rec) < 0) {
      bam_destroy1(rec);
      rec = 0; // Exit after this iteration.
    }

    if((!rec) || (prevrec->core.l_qname && strcmp(bam_get_qname(prevrec), bam_get_qname(rec)))) {

      // Start of a new block. Record the block just concluded.

      if(rec && strlen(bam_get_qname(prevrec)) && strnum_cmp(bam_get_qname(rec), bam_get_qname(prevrec)) < 0) {
	fprintf(stderr, "Input went backwards from %s to %s\n", bam_get_qname(rec), bam_get_qname(prevrec));
	exit(1);
      }

      matchingRecs.sort();

      int blockLim = matchingRecs.recs.size();
      int secondMateBegins = matchingRecs.recs.size();

      for(int i = 0; i < blockLim; ++i) {
	if(flag2mate(matchingRecs.recs[i]) == 2) {
	  secondMateBegins = i;
	  break;
	}
      }

      for(int mate = 1; mate <= 2; ++mate) {

	int startRec = mate == 1 ? 0 : secondMateBegins;
	int limRec = mate == 1 ? secondMateBegins : blockLim;

	if(startRec == limRec)
	  continue;

	contigs.data.clear();
	contigs.data.reserve(limRec - startRec);

	for(int i = startRec; i != limRec; ++i)
	  contigs.data.push_back(matchingRecs.recs[i]->core.tid);

	std::sort(contigs.data.begin(), contigs.data.end());
	std::vector<int>::iterator uniqit = std::unique(contigs.data.begin(), contigs.data.end());
	contigs.data.resize(std::distance(contigs.data.begin(), uniqit));

	++(counts[contigs]);

      }

      matchingRecs.clear();

    }

    if(rec) {
      matchingRecs.copy_add(rec);
      bam_copy1(prevrec, rec);
    }

  } while(rec);

  for(std::map<CmpIntVector, int>::const_iterator it = counts.begin(), itend = counts.end(); it != itend; ++it) {

    fprintf(fo, "%d", it->second);
    for(std::vector<int>::const_iterator contigsit = it->first.data.begin(), contigsend = it->first.data.end();
	contigsit != contigsend; ++contigsit) {

      fprintf(fo, ",%d", *contigsit);

    }

    fprintf(fo, "\n");

  }

  hts_close(hfi);
  fclose(fo);

}
