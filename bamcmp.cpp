
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <vector>
#include <algorithm>
#include <string>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

enum cmptypes {
  
  cmptype_sequence,
  cmptype_mate

};

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

static void usage() {

  fprintf(stderr, "Usage: intersect -a input1.s/b/cram -b input2.s/b/cram [-m matching_records.s/b/cram] [-1 first_only.s/b/cram] [-2 second_only.s/b/cram] [-t nthreads] [-B iobuffersize] [-c comparison_type]\n");
  exit(1);

}

static htsFile* hts_begin_or_die(const char* filename, const char* mode, bam_hdr_t* header, int nthreads, size_t buffersize) {

  htsFile* hf = hts_open(filename, mode);
  if(!hf) {
    fprintf(stderr, "Failed to open %s\n", filename);
    exit(1);
  }

  // header should be 0 if we're opening to read
  if(header && sam_hdr_write(hf, header)) {
    fprintf(stderr, "Failed to write header for %s\n", filename);
    exit(1);
  }

  if(nthreads != 1) {
    hts_set_threads(hf, nthreads);
    if(hf->is_bin && strchr(mode, 'r')) {
      // Cache space required for BAM readahead
      bgzf_set_cache_size(hf->fp.bgzf, BGZF_MAX_BLOCK_SIZE * nthreads * 256);
    }
  }
  if(buffersize)
    hts_set_opt(hf, HTSOL_FILEIO, HTS_FILEIO_BUFFER_SIZE, buffersize);

  return hf;

}

struct SamReader {

  htsFile* hf;
  bam_hdr_t* header;
  bam1_t *rec;
  bool eof;
  SamReader(htsFile* _hf, bam_hdr_t* _header) : hf(_hf), header(_header), eof(false) {
    rec = bam_init1();
    next();
  }

  ~SamReader() {
    bam_destroy1(rec);
  }

  bool is_eof() const {
    return eof;
  }

  void next() {
    if(eof)
      return;
    if(sam_read1(hf, header, rec) < 0)
      eof = true;
  }

};

cmptypes cmptype;

static int32_t get_alignment_score(bam1_t* rec) {
  
  uint8_t* score_rec = bam_aux_get(rec, "AS");
  if(!score_rec) {
    return rec->core.qual;
  }
  else {
    return bam_aux2i(score_rec);
  }

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

static bool sequence_eq(const bam1_t* a, const bam1_t* b) {

  const unsigned char* seqa = bam_get_seq(a);
  const unsigned char* seqb = bam_get_seq(b);
  if(a->core.l_qseq != b->core.l_qseq)
    return false;

  for(unsigned i = 0, ilim = a->core.l_qseq; i != ilim; ++i) 
    if(bam_seqi(seqa, i) != bam_seqi(seqb, i))
      return false;

  return true;
    
}

static bool bamrec_eq(const bam1_t* a, const bam1_t* b) {
  switch(cmptype) {
  case cmptype_sequence:
    return sequence_eq(a, b);
  case cmptype_mate:
    return mate_eq(a, b);
  }
}

static bool mate_lt(const bam1_t* a, const bam1_t* b) {

  return flag2mate(a) < flag2mate(b);

}

static bool sequence_lt(const bam1_t* a, const bam1_t* b) {

  const unsigned char* seqa = bam_get_seq(a);
  const unsigned char* seqb = bam_get_seq(b);
  for(unsigned i = 0, ilim = std::min(a->core.l_qseq, b->core.l_qseq); i != ilim; ++i) {
    char basea = bam_seqi(seqa, i);
    char baseb = bam_seqi(seqb, i);
    if(basea < baseb)
      return true;
    else if(basea > baseb)
      return false;
  }

  return a->core.l_qseq < b->core.l_qseq;

}

static bool bamrec_lt(const bam1_t* a, const bam1_t* b) {
  switch(cmptype) {
  case cmptype_sequence:
    return sequence_lt(a, b);
  case cmptype_mate:
    return mate_lt(a, b);
  }
}

struct BamRecVector {

  std::vector<bam1_t*> recs;

  BamRecVector() {}
  ~BamRecVector() {
    clear();
  }

  void take_add(bam1_t* src) {
    recs.push_back(src);
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
    std::sort(recs.begin(), recs.end(), bamrec_lt);
  }

};

int main(int argc, char** argv) {

  char *in1_name = 0, *in2_name = 0, *both_name = 0, *first_name = 0, *second_name = 0;

  int nthreads = 1;
  size_t buffersize = 0;
  char* cmptypestr = "sequence";

  char c;
  while ((c = getopt(argc, argv, "a:b:m:1:2:c:t:B:")) >= 0) {
    switch (c) {
    case 'a':
      in1_name = optarg;
      break;
    case 'b':
      in2_name = optarg;
      break;
    case 'm':
      both_name = optarg;
      break;
    case '1':
      first_name = optarg;
      break;
    case '2':
      second_name = optarg;
      break;
    case 't':
      nthreads = atoi(optarg);
      break;
    case 'B':
      buffersize = atoi(optarg);
      break;
    case 'c':
      cmptypestr = optarg;
      break;
    default:
      usage();
    }
  }

  if(!in1_name)
    usage();
  if(!in2_name)
    usage();
  if(!(both_name || first_name || second_name)) {
    fprintf(stderr, "intersect is useless without at least one of -m -1 and -2\n");
    usage();
  }

  if(!strcmp(cmptypestr, "sequence"))
    cmptype = cmptype_sequence;
  else if(!strcmp(cmptypestr, "mate"))
    cmptype = cmptype_mate;
  else
    usage();
  
  htsFile *in1hf = 0, *in2hf = 0, *both_out = 0, *first_out = 0, *second_out = 0;
  in1hf = hts_begin_or_die(in1_name, "r", 0, nthreads, buffersize);
  in2hf = hts_begin_or_die(in2_name, "r", 0, nthreads, buffersize);

  bam_hdr_t* header1 = sam_hdr_read(in1hf);
  bam_hdr_t* header2 = sam_hdr_read(in2hf);

  // TODO: what if the headers disagree?

  if(both_name)
    both_out = hts_begin_or_die(both_name, "wb0", header1, nthreads, buffersize);
  if(first_name)
    first_out = hts_begin_or_die(first_name, "wb0", header1, nthreads, buffersize);
  if(second_name)
    second_out = hts_begin_or_die(second_name, "wb0", header1, nthreads, buffersize);

  SamReader in1(in1hf, header1);
  SamReader in2(in2hf, header2);

  BamRecVector seqs1, seqs2;

  while((!in1.is_eof()) && (!in2.is_eof())) {

    std::string qname1 = bam_get_qname(in1.rec);
    std::string qname2 = bam_get_qname(in2.rec);

    if(qname1 == qname2) {
      
      seqs1.clear();
      seqs2.clear();

      std::string qn;
      while((!in1.is_eof()) && (qn = bam_get_qname(in1.rec)) == qname1) {
	seqs1.copy_add(in1.rec);
	in1.next();
      }
      
      seqs1.sort();

      while((!in2.is_eof()) && (qn = bam_get_qname(in2.rec)) == qname2) {
	seqs2.copy_add(in2.rec);
	in2.next();
      }

      seqs2.sort();

      int idx1 = 0, idx2 = 0;
      while(idx1 < seqs1.recs.size() && idx2 < seqs2.recs.size()) {

	if(bamrec_eq(seqs1.recs[idx1], seqs2.recs[idx2])) {

	  int32_t score1 = get_alignment_score(seqs1.recs[idx1]);
	  int32_t score2 = get_alignment_score(seqs2.recs[idx2]);

	  // Specifically the second file might have multiple mappings. Use the best alignment score given,
	  // and skip the rest.
	  while(idx2 + 1 < seqs2.recs.size() && bamrec_eq(seqs1.recs[idx1], seqs2.recs[idx2 + 1])) {
	    ++idx2;
	    score2 = std::max(score2, get_alignment_score(seqs2.recs[idx2]));
	  }
	  	  
	  uint32_t human_better = 1;
	  
	  if(score1 > score2)
	    bam_aux_append(seqs1.recs[idx1], "hb", 'i', sizeof(uint32_t), (uint8_t*)&human_better);  

	  if(both_out)
	    sam_write1(both_out, header1, seqs1.recs[idx1]);

	  ++idx1; ++idx2;

	}
	else if(bamrec_lt(seqs1.recs[idx1], seqs2.recs[idx2])) {
	  if(first_out)
	    sam_write1(first_out, header1, seqs1.recs[idx1]);
	  ++idx1;
	}
	else {
	  if(second_out)
	    sam_write1(second_out, header1, seqs2.recs[idx2]);
	  ++idx2;
	}
		   
      }

      for(;idx1 < seqs1.recs.size(); ++idx1) 
	if(first_out)
	  sam_write1(first_out, header1, seqs1.recs[idx1]);

      for(;idx2 < seqs2.recs.size(); ++idx2) 
	if(second_out)
	  sam_write1(second_out, header1, seqs2.recs[idx2]);

    }
    else if(strnum_cmp(qname1.c_str(), qname2.c_str()) < 0) {

      if(first_out)
	sam_write1(first_out, header1, in1.rec);
      in1.next();

    }
    else {

      if(second_out)
	sam_write1(second_out, header1, in2.rec);
      in2.next();

    }

  }

  // One or other file has reached EOF. Write the remainder as first- or second-only records.

  if(first_out) {
    while(!in1.is_eof()) {
      sam_write1(first_out, header1, in1.rec);
      in1.next();
    }
  }

  if(second_out) {
    while(!in2.is_eof()) {
      sam_write1(second_out, header1, in2.rec);
      in2.next();
    }
  }

  hts_close(in1hf);
  hts_close(in2hf);
  if(both_out)
    hts_close(both_out);
  if(first_out)
    hts_close(first_out);
  if(second_out)
    hts_close(second_out);

}
