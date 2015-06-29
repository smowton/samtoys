
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

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

enum scoringmethods {
  
  scoringmethod_nmatches,
  scoringmethod_astag,
  scoringmethod_mapq

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

static bool mixed_ordering = true;

static int qname_cmp(const char* qa, const char* qb) {

  if(mixed_ordering)
    return strnum_cmp(qa, qb);
  else
    return strcmp(qa, qb);
  
}

static void usage() {

  fprintf(stderr, "Usage: intersect -a input1.s/b/cram -b input2.s/b/cram [-1 first_only.xam] [-A first_better.xam] [-B second_better.xam] [-2 second_only.xam] [-t nthreads] [-n | -N] [-s scoring_method] [-c comparison_type]\n");
  fprintf(stderr, "\t-n\tExpect input sorted as per samtools -n (Mixed string / integer ordering, default)\n");
  fprintf(stderr, "\t-N\tExpect input sorted as per Picard / htsjdk name ordering (lexical ordering)\n");
  fprintf(stderr, "\t-s match\tScore hits by the number of bases that match the reference, as given by the CIGAR string and NM / MD attributes\n");
  fprintf(stderr, "\t-s as\tScore hits according to the AS attribute written by some aligners\n");
  fprintf(stderr, "\t-s mapq\tScore hits according to the MAPQ SAM field\n");
  exit(1);

}

static htsFile* hts_begin_or_die(const char* filename, const char* mode, bam_hdr_t* header, int nthreads) {

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

  if(nthreads != 1)
    hts_set_threads(hf, nthreads);

  return hf;

}

cmptypes cmptype;
scoringmethods scoringmethod;

bool warned_nm_anomaly = false;
bool warned_nm_md_tags = false;

static bool aux_is_int(uint8_t* rec) {

  switch(*rec) {
  case 'c':
  case 'C':
  case 's':
  case 'S':
  case 'i':
  case 'I':
    return true;
  default:
    return false;
  }

}

static uint32_t get_alignment_score(bam1_t* rec) {
  
  switch(scoringmethod) {
  case scoringmethod_nmatches:
    {

      bool seen_equal_or_diff = false;
      int32_t cigar_total = 0;
      const uint32_t* cigar = bam_get_cigar(rec);

      int32_t indel_edit_distance = 0;

      for(int i = 0; i < rec->core.n_cigar; ++i) {

	// CIGAR scoring: score points for matching bases, and negatives for deletions
	// since otherwise 10M10D10M would score the same as 20M. Insertions, clipping etc
	// don't need to score a penalty since they skip bases in the query.
	// CREF_SKIP (N / intron-skip operator) is acceptable: 10M1000N10M is as good as 20M.
	// Insertions are counted to correct the NM tag below only.
	
	int32_t n = bam_cigar_oplen(cigar[i]);
	switch(bam_cigar_op(cigar[i])) {

	case BAM_CEQUAL:
	  seen_equal_or_diff = true;
	  // fall through
	case BAM_CMATCH:
	  cigar_total += n;
	  break;

	case BAM_CDEL:
	  indel_edit_distance += n;
	  cigar_total -= n;
	  break;

	case BAM_CDIFF:
	  seen_equal_or_diff = true;
	  break;

	case BAM_CINS:
	  indel_edit_distance += n;
	  break;

	default:
	  break;

	}

      }

      // The BAM_CMATCH operator (unlike BAM_CEQUAL or BAM_CDIFF) could mean a match or a mismatch
      // with same length (e.g. a SNP). If the file doesn't seem to use the advanced operators try to
      // spot mismatches from metadata tags.

      if(!seen_equal_or_diff) {
	uint8_t* nm_rec = bam_aux_get(rec, "NM");
	if(nm_rec && aux_is_int(nm_rec)) {
	  int32_t nm = bam_aux2i(nm_rec);
	  if(nm < indel_edit_distance) {
	    if(!warned_nm_anomaly) {
	      fprintf(stderr, "Warning: anomaly in record %s: NM is %d but there are at least %d indel bases in the CIGAR string\n", bam_get_qname(rec), nm, indel_edit_distance);
	      fprintf(stderr, "There may be more records with this problem, but the warning will not be repeated\n");
	      warned_nm_anomaly = true;
	    }
	  }
	  else {
	    cigar_total -= (bam_aux2i(nm_rec) - indel_edit_distance);
	    seen_equal_or_diff = true;
	  }
	}
      }

      if(!seen_equal_or_diff) {
	uint8_t* md_rec = bam_aux_get(rec, "MD");
	if(md_rec) {

	  char* mdstr = bam_aux2Z(md_rec);
	  if(mdstr) {
	    
	    seen_equal_or_diff = true;

	    bool in_deletion = false;
	    for(; *mdstr; ++mdstr) {
	      
	      // Skip deletions, which are already penalised.
	      // Syntax seems to be: numbers mean base strings that match the reference; ^ followed by letters means
	      // a deletion; letters without the preceding ^ indicate a mismatch.

	      char c = *mdstr;
	      if(c == '^')
		in_deletion = true;
	      else if(isdigit(c))
		in_deletion = false;
	      else if(!in_deletion) {
		// Mismatch
		cigar_total--;
	      }

	    }

	  }

	}
      }

      if((!seen_equal_or_diff) && !warned_nm_md_tags) {

	fprintf(stderr, "Warning: input file does not use the =/X CIGAR operators, or include NM or MD tags, so I have no way to spot length-preserving reference mismatches.\n");
	fprintf(stderr, "At least record %s exhibited this problem; there may be others but the warning will not be repeated. I will assume M CIGAR operators indicate a match.\n", bam_get_qname(rec));
	warned_nm_md_tags = true;

      }

      return std::max(cigar_total, 0);

    }
    // End the CIGAR string scoring method. Thankfully the others are much simpler to implement:

  case scoringmethod_astag:
    {
      uint8_t* score_rec = bam_aux_get(rec, "AS");
      if(!score_rec) {
	fprintf(stderr, "Fatal: At least record %s doesn't have an AS tag as required.\n", bam_get_qname(rec));
	exit(1);
      }
      return bam_aux2i(score_rec);
    }

  case scoringmethod_mapq:
    return rec->core.qual;
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

struct SamReader {

  htsFile* hf;
  bam_hdr_t* header;
  bam1_t *prev_rec;
  bam1_t *rec;
  bool eof;
  std::string filename;

  SamReader(htsFile* _hf, bam_hdr_t* _header, const char* fname) : hf(_hf), header(_header), eof(false), filename(fname) {
    rec = bam_init1();
    prev_rec = bam_init1();
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

    if(rec->data)
      bam_copy1(prev_rec, rec);

    if(sam_read1(hf, header, rec) < 0)
      eof = true;
    
    if(prev_rec->data && (!eof) && qname_cmp(bam_get_qname(rec), bam_get_qname(prev_rec)) < 0) {
      fprintf(stderr, "Order went backwards! In file %s, record %s belongs before %s. Re-sort your files and try again.\n", filename.c_str(), bam_get_qname(rec), bam_get_qname(prev_rec));
      if(mixed_ordering)
	fprintf(stderr, "Expected order was the mixed string/integer ordering produced by samtools sort -n; use -N to switch to Picard / htsjdk string ordering\n");
      else
	fprintf(stderr, "Expected order was Picard / htsjdk string ordering; use -n to switch to samtools sort -n ordering\n");
      exit(1);
    }
      
  }

};

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

  char *in1_name = 0, *in2_name = 0, *firstbetter_name = 0, *secondbetter_name = 0, *first_name = 0, *second_name = 0;

  int nthreads = 1;
  size_t buffersize = 0;
  const char* cmptypestr = "sequence";
  const char* scoring_method_string = "match";

  char c;
  while ((c = getopt(argc, argv, "a:b:m:1:2:c:t:A:B:nNs:")) >= 0) {
    switch (c) {
    case 'a':
      in1_name = optarg;
      break;
    case 'b':
      in2_name = optarg;
      break;
    case '1':
      first_name = optarg;
      break;
    case '2':
      second_name = optarg;
      break;
    case 'A':
      firstbetter_name = optarg;
      break;
    case 'B':
      secondbetter_name = optarg;
      break;
    case 't':
      nthreads = atoi(optarg);
      break;
    case 'c':
      cmptypestr = optarg;
      break;
    case 'n':
      mixed_ordering = true;
      break;
    case 'N':
      mixed_ordering = false;
      break;
    case 's':
      scoring_method_string = optarg;
      break;
    default:
      usage();
    }
  }

  if(!in1_name)
    usage();
  if(!in2_name)
    usage();
  if(!(first_name || second_name || firstbetter_name || secondbetter_name)) {
    fprintf(stderr, "intersect is useless without at least one of -1, -2, -A or -B\n");
    usage();
  }

  if(!strcmp(scoring_method_string, "match"))
    scoringmethod = scoringmethod_nmatches;
  else if(!strcmp(scoring_method_string, "mapq"))
    scoringmethod = scoringmethod_mapq;
  else if(!strcmp(scoring_method_string, "as"))
    scoringmethod = scoringmethod_astag;
  else
    usage();

  if(!strcmp(cmptypestr, "sequence"))
    cmptype = cmptype_sequence;
  else if(!strcmp(cmptypestr, "mate"))
    cmptype = cmptype_mate;
  else
    usage();
  
  htsFile *in1hf = 0, *in2hf = 0, *firstbetter_out = 0, *secondbetter_out = 0, *first_out = 0, *second_out = 0;
  in1hf = hts_begin_or_die(in1_name, "r", 0, nthreads);
  in2hf = hts_begin_or_die(in2_name, "r", 0, nthreads);

  bam_hdr_t* header1 = sam_hdr_read(in1hf);
  bam_hdr_t* header2 = sam_hdr_read(in2hf);

  // Permit the outputs using like headers to share a file if they gave the same name.
  
  if(firstbetter_name)
    firstbetter_out = hts_begin_or_die(firstbetter_name, "wb0", header1, nthreads);
  if(secondbetter_name)
    secondbetter_out = hts_begin_or_die(secondbetter_name, "wb0", header2, nthreads);
  if(first_name) {
    if(firstbetter_name && !strcmp(firstbetter_name, first_name))
      first_out = firstbetter_out;
    else
      first_out = hts_begin_or_die(first_name, "wb0", header1, nthreads);
  }
  if(second_name) {
    if(secondbetter_name && !strcmp(secondbetter_name, second_name))
      second_out = secondbetter_out;
    else
      second_out = hts_begin_or_die(second_name, "wb0", header2, nthreads);
  }

  SamReader in1(in1hf, header1, in1_name);
  SamReader in2(in2hf, header2, in2_name);

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

	  uint32_t score1 = get_alignment_score(seqs1.recs[idx1]);
	  uint32_t score2 = get_alignment_score(seqs2.recs[idx2]);

	  // Either input may have multiple candidate matches. Compare the best match found in each group
	  // and then emit the whole group as firstbetter or secondbetter.
	 
	  int group_start_idx1 = idx1, group_start_idx2 = idx2;

	  while(idx1 + 1 < seqs1.recs.size() && bamrec_eq(seqs1.recs[group_start_idx1], seqs1.recs[idx1 + 1])) {
	    ++idx1;
	    score1 = std::max(score1, get_alignment_score(seqs1.recs[idx1]));
	  }

	  while(idx2 + 1 < seqs2.recs.size() && bamrec_eq(seqs1.recs[group_start_idx1], seqs2.recs[idx2 + 1])) {
	    ++idx2;
	    score2 = std::max(score2, get_alignment_score(seqs2.recs[idx2]));
	  }
	  	  
	  BamRecVector* writeseqs;
	  uint32_t writefrom, writeto;
	  htsFile* writefile;
	  bam_hdr_t* writeheader;

	  if(score1 > score2) {
	    writeseqs = &seqs1;
	    writefrom = group_start_idx1;
	    writeto = idx1;
	    writefile = firstbetter_out;
	    writeheader = header1;
	  }
	  else {
	    writeseqs = &seqs2;
	    writefrom = group_start_idx2;
	    writeto = idx2;
	    writefile = secondbetter_out;
	    writeheader = header2;
	  }
	 
	  if(writefile) {

	    for(; writefrom <= writeto; ++writefrom) {
	      
	      bam_aux_append(writeseqs->recs[writefrom], "as", 'i', sizeof(uint32_t), (uint8_t*)&score1);  
	      bam_aux_append(writeseqs->recs[writefrom], "bs", 'i', sizeof(uint32_t), (uint8_t*)&score2);  
	      sam_write1(writefile, writeheader, writeseqs->recs[writefrom]);
	      
	    }

	  }

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
    else if(qname_cmp(qname1.c_str(), qname2.c_str()) < 0) {

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
  if(first_out)
    hts_close(first_out);
  if(second_out)
    hts_close(second_out);

  // Only close these if the hts handle wasn't shared with first/second_out.
  if(firstbetter_out && ((!first_out) || strcmp(first_name, firstbetter_name) != 0))
    hts_close(firstbetter_out);
  if(secondbetter_out && ((!second_out) || strcmp(second_name, secondbetter_name) != 0))
    hts_close(secondbetter_out);

}
