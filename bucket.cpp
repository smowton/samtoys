
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/ksort.h>

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

// Stolen from bam_sort.c to sort by qname.

typedef bam1_t *bam1_p;

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

static inline int bam1_lt(const bam1_p a, const bam1_p b)
{
  int t = strnum_cmp(bam_get_qname(a), bam_get_qname(b));
  return (t < 0 || (t == 0 && (a->core.flag&0xc0) < (b->core.flag&0xc0)));
}
KSORT_INIT(sort, bam1_p, bam1_lt)

int main(int argc, char *argv[])
{

  const char* tmpprefix = argv[1];
	
  samFile* in = sam_open("-", "r");
  if(!in) {
    fprintf(stderr, "Open\n");
    exit(1);
  }

  bam_hdr_t* header = sam_hdr_read(in);

  int n_outs = argc - 1;
  samFile** outs = (samFile**)calloc(sizeof(samFile*), n_outs);
  if(!outs) {
    fprintf(stderr, "Malloc 1\n");
    exit(1);
  }
	
  char buf[4096];
  int i;

  for(i = 0; i < n_outs; ++i) {
    snprintf(buf, 4096, "%s.%d.bam", tmpprefix, i);
    buf[4095] = '\0';
    outs[i] = sam_open(buf, "wb0");
    if(!outs[i]) {
      fprintf(stderr, "Open out %d\n", i);
      exit(1);
    }
    sam_hdr_write(outs[i], header);
  }

  bam1_t* buckets = (bam1_t*)calloc(sizeof(bam1_t), n_outs - 1);
  if(!buckets) {
    fprintf(stderr, "Malloc 2\n");
    exit(1);
  }

  for(i = 0; i < (n_outs - 1); ++i) {

    char contig[128];
    int idx = 0;

    // Parse bucket start points specified as contigname:index
    int r = sscanf(argv[i + 2], "%127[a-zA-Z0-9]:%d", contig, &idx);
    contig[127] = '\0';

    if(r != 2 || !idx) {
      fprintf(stderr, "Malformed argument %s\n", argv[i + 2]);
      exit(1);
    }

    // Only these fields are relevant to bam entry comparison
    buckets[i].core.tid = bam_name2id(header, contig);
    buckets[i].core.pos = idx;
    // Set on forward strand
    buckets[i].core.flag = 0;

    if(buckets[i].core.tid == -1) {
      fprintf(stderr, "Bad contig %s\n", contig);
      exit(1);
    }

  }

  bam1_t rd;
	
  while(sam_read1(in, header, &rd) >= 0) {

    // I expect a small number of buckets. Switch to binary search if larger
    // scale usage is desired.

    int bucket = 0;
    while(bucket < (n_outs - 1) && bam1_lt(&buckets[bucket], &rd))
      ++bucket;

    sam_write1(outs[bucket], header, &rd);

  }

  sam_close(in);
  for(i = 0; i < n_outs; ++i)
    sam_close(outs[i]);

  return 0;

}
