

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <vector>
#include <algorithm>

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

static int get_as(const bam1_t* rec) {

  uint8_t* score_rec = bam_aux_get(rec, "AS");
  if(!score_rec) {
    fprintf(stderr, "Fatal: At least record %s doesn't have an AS tag as required.\n", bam_get_qname(rec));
    exit(1);
  }
  return bam_aux2i(score_rec);

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

int main(int argc, char** argv) {
  
  if(argc < 4 || argc > 5) {

    fprintf(stderr, "Usage: filter_hits in.xam out.xam maxhits [trim]\n");
    exit(1);

  }

  htsFile* hfi = hts_open(argv[1], "r");
  if(!hfi) {
    fprintf(stderr, "Failed to open %s\n", argv[1]);
    exit(1);
  }

  htsFile* hfo = hts_open(argv[2], "wb");
  if(!hfo) {
    fprintf(stderr, "Failed to open %s\n", argv[2]);
    exit(1);
  }

  int maxhits = atoi(argv[3]);
  if(maxhits <= 0) {
    fprintf(stderr, "Argument 3 must be a number >= 1 specifying how many hits we're willing to keep\n");
    exit(1);
  }

  bool do_trim = false;
  if(argc == 5) {
    if(!strcmp(argv[4], "trim"))
      do_trim = true;
    else {
      fprintf(stderr, "Argument 4 should be 'trim' if present"); 
      exit(1);
    }
  }

  bam_hdr_t* header = sam_hdr_read(hfi);
  sam_hdr_write(hfo, header);

  bam1_t* prevrec = bam_init1();
  bam1_t* rec = bam_init1();

  std::vector<int64_t> counts;

  BamRecVector matchingRecs;

  do {

    if(sam_read1(hfi, header, rec) < 0) {
      bam_destroy1(rec);
      rec = 0; // Exit after this iteration.
    }

    if((!rec) || (prevrec->core.l_qname && strcmp(bam_get_qname(prevrec), bam_get_qname(rec)))) {

      // Start of a new block. Write the records in the previous block out if it's small enough.

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

	int blockSize = limRec - startRec;

	if(blockSize >= counts.size())
	  counts.resize(blockSize + 1);
	++(counts[blockSize]);

	if(do_trim) {
	  
	  if(maxhits != 1) {
	    fprintf(stderr, "Trimming with maxhits != 1 not implemented yet\n");
	    exit(1);
	  }

	  int bestRec = -1;
	  int bestAS = -1;

	  for(int i = startRec; i != limRec; ++i) {
	   
	    int as = get_as(matchingRecs.recs[i]);
	    if(bestAS < as) {
	      bestAS = as;
	      bestRec = i;
	    }

	  }

	  if(bestRec != -1) {

	    startRec = bestRec;
	    limRec = startRec + maxhits;
	    blockSize = maxhits;

	  }

	}

	if(blockSize <= maxhits) {

	  for(int i = startRec; i != limRec; ++i)
	    sam_write1(hfo, header, matchingRecs.recs[i]);

	}

      }

      matchingRecs.clear();

    }

    if(rec) {
      matchingRecs.copy_add(rec);
      bam_copy1(prevrec, rec);
    }

  } while(rec);

  fprintf(stderr, "Counts histogram:\n");

  for(int i = 1; i < counts.size(); ++i)
    fprintf(stderr, "%d: %ld\n", i, counts[i]);

  hts_close(hfi);
  hts_close(hfo);

}
