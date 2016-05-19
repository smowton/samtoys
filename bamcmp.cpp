
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

enum scoringmethods {
  
  scoringmethod_nmatches,
  scoringmethod_astag,
  scoringmethod_mapq,
  scoringmethod_balwayswins

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

  fprintf(stderr, "Usage: intersect -1 input1.s/b/cram -2 input2.s/b/cram [-a first_only.xam] [-b second_only.xam] [-A first_better.xam] [-B second_better.xam] [-C first_worse.xam] [-D second_worse.xam] [-t nthreads] [-n | -N] [-s scoring_method]\n");
  fprintf(stderr, "\t-n\tExpect input sorted as per samtools -n (Mixed string / integer ordering, default)\n");
  fprintf(stderr, "\t-N\tExpect input sorted as per Picard / htsjdk name ordering (lexical ordering)\n");
  fprintf(stderr, "\t-s match\tScore hits by the number of bases that match the reference, as given by the CIGAR string and NM / MD attributes (default)\n");
  fprintf(stderr, "\t-s as\tScore hits according to the AS attribute written by some aligners\n");
  fprintf(stderr, "\t-s mapq\tScore hits according to the MAPQ SAM field\n");
  fprintf(stderr, "\t-s balwayswins\tAlways award hits to input B, regardless of alignment scores (equivalent to filtering A by any read mapped in B)\n");
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

class htsFileWrapper {

  std::string fname;
  const char* mode;
  htsFile* hts;
  uint32_t refCount;
  int nthreads;
  int header2_offset;

  bam_hdr_t* header1;
  bam_hdr_t* header2;
  bam_hdr_t* headerout;

  void checkHeaderNotWritten() {
    if(headerout) {
      fprintf(stderr, "htsFileWrapper header changed after already being written\n");
      exit(1);
    }
  }

public:
  
  void setHeader1(bam_hdr_t* h1) {
    checkHeaderNotWritten();
    header1 = h1;
  }

  void setHeader2(bam_hdr_t* h2) {
    checkHeaderNotWritten();
    header2 = h2;
  }

  void ref() {
    ++refCount;
  }

  uint32_t unref() {

    checkStarted();

    if(refCount == 1)
      hts_close(hts);

    return --refCount;

  }

  htsFileWrapper(const std::string& _fname, const char* _mode, int _nthreads) : 
    fname(_fname), mode(_mode), hts(0), refCount(1), nthreads(_nthreads), header2_offset(0), header1(0), header2(0), headerout(0) { }

  void checkStarted() {

    if(hts)
      return;
    
    if(!(header1 || header2)) {
      fprintf(stderr, "Started writing records without any header\n");
      exit(1);
    }
    else if(!header2)
      headerout = header1;
    else if(!header1)
      headerout = header2;
    else {

      header2_offset = header1->n_targets;

      headerout = bam_hdr_dup(header1);
      headerout->n_targets += header2->n_targets;
      headerout->target_len = (uint32_t*)realloc(headerout->target_len, sizeof(uint32_t) * headerout->n_targets);
      headerout->target_name = (char**)realloc(headerout->target_name, sizeof(char*) * headerout->n_targets);

      if((!headerout->target_len) || (!headerout->target_name))
	goto oom;

      // Prefix header1 names with A_

      for(int i = 0; i < header1->n_targets; ++i) {
	int len = strlen(headerout->target_name[i]);
	headerout->target_name[i] = (char*)realloc(headerout->target_name[i], len + 3); // A_ prefix plus null terminator
	if(!headerout->target_name[i])
	  goto oom;
	memmove(headerout->target_name[i] + 2, headerout->target_name[i], len + 1);
	memcpy(headerout->target_name[i], "A_", 2);
      }
      
      // Copy header2 names; add B_ prefix

      for(int i = 0; i < header2->n_targets; ++i) {
	headerout->target_len[i + header2_offset] = header2->target_len[i];
	headerout->target_name[i + header2_offset] = (char*)malloc(strlen(header2->target_name[i]) + 3);
	if(!headerout->target_name[i + header2_offset])
	  goto oom;
	sprintf(headerout->target_name[i + header2_offset], "B_%s", header2->target_name[i]);
      }

      // Find all text @SQ records in header2:
      
      std::vector<std::pair<int, int> > sqheaders;

      int offset = 0;
      char* sq_start;
      while((sq_start = strstr(header2->text + offset, "@SQ")) != 0) {

	// Must occur at the start of the header or after a newline:
	if(sq_start != header2->text && (*(sq_start - 1)) != '\n') {
	  ++offset;
	  continue;
	}

	int start_off = sq_start - header2->text;
	int end_off;

	char* nl = strchr(sq_start, '\n');
	if(!nl)
	  end_off = strlen(header2->text);
	else
	  end_off = nl - header2->text;

	sqheaders.push_back(std::make_pair(start_off, end_off));
	
	offset = end_off;

      }

      // Count how much extra we're going to add to header1:

      int additional_bytes = (2 * header1->n_targets); // A_ prefix for all existing @SQs
      
      for(int i = 0, ilim = sqheaders.size(); i != ilim; ++i)
	additional_bytes += ((sqheaders[i].second - sqheaders[i].first) + 3); // Existing header, plus a newline, plus B_ prefix.

      free(headerout->text);
      headerout->text = (char*)malloc(headerout->l_text + additional_bytes + 1); // Null terminator as well?
      headerout->l_text += additional_bytes;
      if(!headerout->text)
	goto oom;

      // Rewrite: copy header1 lines, except for @SQs which gain an A_ prefix, and insert the (prefixed) header2 @SQs
      // after the last one from header1.

      char* nl;
      char* read_offset = header1->text;
      char* write_offset = headerout->text;
      bool more_to_read = true;
      int sqs_remaining = header1->n_targets;
      while(more_to_read) {

	nl = strchr(read_offset, '\n');
	if(!nl) {
	  nl = header1->text + header1->l_text;
	  more_to_read = false;
	}
	else {
	  // Include the newline when copying.
	  ++nl;
	}

	if(strncmp(read_offset, "@SQ", 3) == 0) {

	  char* sn_offset = strstr(read_offset, "SN:");
	  if(sn_offset) {
	    
	    // Prefix A_ to the @SQ line
	    char* insert_offset = sn_offset + 3;
	    int len = insert_offset - read_offset;
	    memcpy(write_offset, read_offset, len);

	    write_offset += len; read_offset += len;

	    memcpy(write_offset, "A_", 2);
	    write_offset += 2;

	    len = nl - read_offset;
	    memcpy(write_offset, read_offset, len);

	    write_offset += len; read_offset += len;

	    if(!more_to_read) {
	      // Adding at the end; add newline.
	      *(write_offset++) = '\n';
	    }

	    if(!(--sqs_remaining)) {

	      // That was the last of header1's @SQ lines. Insert header2's lines here.
	      for(int i = 0, ilim = sqheaders.size(); i != ilim; ++i) {

		int h2readoff = sqheaders[i].first;
		int h2readlim = sqheaders[i].second;

		char* h2sn = strstr(header2->text + h2readoff, "SN:");
		if(h2sn) {

		  // Skip over the SN: tag
		  h2sn += 3;

		  int h2insoff = h2sn - header2->text;
		  int len = h2insoff - h2readoff;

		  memcpy(write_offset, header2->text + h2readoff, len);
		  write_offset += len;

		  memcpy(write_offset, "B_", 2);
		  write_offset += 2;

		  len = h2readlim - h2insoff;
		  memcpy(write_offset, header2->text + h2insoff, len);
		  write_offset += len;

		}
		else {

		  int len = h2readlim - h2readoff;
		  memcpy(write_offset, header2->text + h2readoff, len);
		  write_offset += len;

		}

		if(more_to_read || i != ilim - 1)
		  *(write_offset++) = '\n';

	      }

	    }

	    continue;

	  }

	  // Else fall through to just copy the whole line:
	  
	}

	// Just a regular header line: copy it all:
	
	int len = nl - read_offset;
	memcpy(write_offset, read_offset, len);
	write_offset += len; read_offset += len;

      }

      headerout->text[headerout->l_text] = '\0';
				
    }

    // Header complete, now open and write it:
    
    hts = hts_begin_or_die(fname.c_str(), mode, headerout, nthreads);
    return;

  oom:

    fprintf(stderr, "Malloc failure while building combined header\n");
    exit(1);

  }

  void write1(int headerNum, bam1_t* rec) {

    checkStarted();

    if(headerNum == 2) {
      if(rec->core.tid != -1)
	rec->core.tid += header2_offset;
      if(rec->core.mtid != -1)
	rec->core.mtid += header2_offset;
    }

    sam_write1(hts, headerout, rec);

    if(headerNum == 2) {
      if(rec->core.tid != -1)
	rec->core.tid -= header2_offset;
      if(rec->core.mtid != -1)
	rec->core.mtid -= header2_offset;
    }

  }

};

static void htswrapper_close(htsFileWrapper* f) {

  if(!f->unref())
    delete f;

}

static std::vector<std::pair<std::string, htsFileWrapper*> > openOutputs;

static htsFileWrapper* htswrapper_begin_or_die(const char* fname, const char* mode, bam_hdr_t* header, int inputNumber, int nthreads) {

  if(inputNumber != 1 && inputNumber != 2) {
    fprintf(stderr, "inputNumber must be 1 or 2\n");
    exit(1);
  }

  htsFileWrapper* ret = 0;

  std::string sfname(fname);
  for(std::vector<std::pair<std::string, htsFileWrapper*> >::iterator it = openOutputs.begin(), itend = openOutputs.end(); it != itend && !ret; ++it) {

    if(it->first == sfname) {
      it->second->ref();
      ret = it->second;
    }

  }

  if(!ret) {
    ret = new htsFileWrapper(sfname, mode, nthreads);
    openOutputs.push_back(std::make_pair(sfname, ret));
  }

  if(inputNumber == 1)
    ret->setHeader1(header);
  else
    ret->setHeader2(header);  

  return ret;

}

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

static uint32_t get_alignment_score(bam1_t* rec, bool is_input_a) {
  
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

  case scoringmethod_balwayswins:

    // Mapped B records beat any A record, beats an unmapped B record.
    if(is_input_a)
      return 1;
    else if(!(rec->core.flag & BAM_FUNMAP))
      return 2;
    else
      return 0;

  }

}

static int flag2mate(const bam1_t* rec) {

  if(rec->core.flag & BAM_FREAD1)
    return 1;
  else if(rec->core.flag & BAM_FREAD2)
    return 2;
  return 0;

}

static bool bamrec_eq(const bam1_t* a, const bam1_t* b) {
  return flag2mate(a) == flag2mate(b);
}

static bool bamrec_lt(const bam1_t* a, const bam1_t* b) {
  return flag2mate(a) < flag2mate(b);
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

static bool uniqueValue(const std::vector<htsFileWrapper*>& in) {

  bool outValid = false;
  htsFileWrapper* out = 0;

  for(std::vector<htsFileWrapper*>::const_iterator it = in.begin(), itend = in.end(); it != itend; ++it) {

    if(!outValid) {
      out = *it;
      outValid = true;
    }
    else if(out != *it) {
      return false;
    }

  }

  return outValid;

}

static void clearMateInfo(BamRecVector& v) {

  for(int i = 0, ilim = v.recs.size(); i != ilim; ++i) {
    
    uint32_t maten = (uint32_t)flag2mate(v.recs[i]);
    bam_aux_append(v.recs[i], "om", 'i', sizeof(uint32_t), (uint8_t*)&maten);  

    v.recs[i]->core.flag &= ~(BAM_FPROPER_PAIR | BAM_FMREVERSE | BAM_FPAIRED | BAM_FMUNMAP | BAM_FREAD1 | BAM_FREAD2);
    v.recs[i]->core.mtid = -1;
    v.recs[i]->core.mpos = -1;

  }

}

int main(int argc, char** argv) {

  char *in1_name = 0, *in2_name = 0, *firstbetter_name = 0, *secondbetter_name = 0, 
    *firstworse_name = 0, *secondworse_name = 0, *first_name = 0, *second_name = 0;

  int nthreads = 1;
  size_t buffersize = 0;
  const char* cmptypestr = "sequence";
  const char* scoring_method_string = "match";

  char c;
  while ((c = getopt(argc, argv, "a:b:m:1:2:t:A:B:C:D:nNs:")) >= 0) {
    switch (c) {
    case '1':
      in1_name = optarg;
      break;
    case '2':
      in2_name = optarg;
      break;
    case 'a':
      first_name = optarg;
      break;
    case 'b':
      second_name = optarg;
      break;
    case 'A':
      firstbetter_name = optarg;
      break;
    case 'B':
      secondbetter_name = optarg;
      break;
    case 'C':
      firstworse_name = optarg;
      break;
    case 'D':
      secondworse_name = optarg;
      break;
    case 't':
      nthreads = atoi(optarg);
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
  else if(!strcmp(scoring_method_string, "balwayswins"))
    scoringmethod = scoringmethod_balwayswins;
  else
    usage();

  htsFile *in1hf = 0, *in2hf = 0;
  htsFileWrapper *firstbetter_out = 0, *secondbetter_out = 0, *firstworse_out = 0, *secondworse_out = 0, *first_out = 0, *second_out = 0;
  in1hf = hts_begin_or_die(in1_name, "r", 0, nthreads);
  in2hf = hts_begin_or_die(in2_name, "r", 0, nthreads);

  bam_hdr_t* header1 = sam_hdr_read(in1hf);
  bam_hdr_t* header2 = sam_hdr_read(in2hf);

  // Permit the outputs using like headers to share a file if they gave the same name.

  if(firstbetter_name)
    firstbetter_out = htswrapper_begin_or_die(firstbetter_name, "wb0", header1, 1, nthreads);
  if(secondbetter_name)
    secondbetter_out = htswrapper_begin_or_die(secondbetter_name, "wb0", header2, 2, nthreads);
  if(firstworse_name)
    firstworse_out = htswrapper_begin_or_die(firstworse_name, "wb0", header1, 1, nthreads);
  if(secondworse_name)
    secondworse_out = htswrapper_begin_or_die(secondworse_name, "wb0", header2, 2, nthreads);
  if(first_name)
    first_out = htswrapper_begin_or_die(first_name, "wb0", header1, 1, nthreads);
  if(second_name)
    second_out = htswrapper_begin_or_die(second_name, "wb0", header2, 2, nthreads);

  SamReader in1(in1hf, header1, in1_name);
  SamReader in2(in2hf, header2, in2_name);

  BamRecVector seqs1, seqs2;
  std::vector<htsFileWrapper*> seqs1Files, seqs2Files;

  while((!in1.is_eof()) && (!in2.is_eof())) {

    std::string qname1 = bam_get_qname(in1.rec);
    std::string qname2 = bam_get_qname(in2.rec);

    if(qname1 == qname2) {
      
      seqs1.clear();
      seqs2.clear();
      seqs1Files.clear();
      seqs2Files.clear();

      std::string qn;
      while((!in1.is_eof()) && (qn = bam_get_qname(in1.rec)) == qname1) {
	seqs1.copy_add(in1.rec);
	in1.next();
      }
      
      seqs1.sort();
      seqs1Files.resize(seqs1.recs.size(), 0);

      while((!in2.is_eof()) && (qn = bam_get_qname(in2.rec)) == qname2) {
	seqs2.copy_add(in2.rec);
	in2.next();
      }

      seqs2.sort();
      seqs2Files.resize(seqs2.recs.size(), 0);

      int idx1 = 0, idx2 = 0;
      while(idx1 < seqs1.recs.size() && idx2 < seqs2.recs.size()) {

	if(bamrec_eq(seqs1.recs[idx1], seqs2.recs[idx2])) {

	  uint32_t score1 = 0;
	  uint32_t score2 = 0; 

	  int group_start_idx1 = idx1, group_start_idx2 = idx2;

	  score1 = get_alignment_score(seqs1.recs[idx1], true);
	  score2 = get_alignment_score(seqs2.recs[idx2], false);

	  // Either input may have multiple candidate matches. Compare the best match found in each group
	  // and then emit the whole group as firstbetter or secondbetter.
	 
	  while(idx1 + 1 < seqs1.recs.size() && bamrec_eq(seqs1.recs[group_start_idx1], seqs1.recs[idx1 + 1])) {
	    ++idx1;
	    score1 = std::max(score1, get_alignment_score(seqs1.recs[idx1], true));
	  }

	  while(idx2 + 1 < seqs2.recs.size() && bamrec_eq(seqs1.recs[group_start_idx1], seqs2.recs[idx2 + 1])) {
	    ++idx2;
	    score2 = std::max(score2, get_alignment_score(seqs2.recs[idx2], false));
	  }
	  	  
	  for(uint32_t i = group_start_idx1; i <= idx1; ++i) {
	    bam_aux_append(seqs1.recs[i], "as", 'i', sizeof(uint32_t), (uint8_t*)&score1);  
	    bam_aux_append(seqs1.recs[i], "bs", 'i', sizeof(uint32_t), (uint8_t*)&score2);  
	  }

	  for(uint32_t i = group_start_idx2; i <= idx2; ++i) {
	    bam_aux_append(seqs2.recs[i], "as", 'i', sizeof(uint32_t), (uint8_t*)&score1);  
	    bam_aux_append(seqs2.recs[i], "bs", 'i', sizeof(uint32_t), (uint8_t*)&score2);  
	  }

	  htsFileWrapper *firstRecordsFile, *secondRecordsFile;

	  if(score1 > score2) {
	    firstRecordsFile = firstbetter_out;
	    secondRecordsFile = secondworse_out;
	  }
	  else {
	    firstRecordsFile = firstworse_out;
	    secondRecordsFile = secondbetter_out;
	  }
	    
	  for(uint32_t i = group_start_idx1; i <= idx1; ++i)
	    seqs1Files[i] = firstRecordsFile;

	  for(uint32_t i = group_start_idx2; i <= idx2; ++i)
	    seqs2Files[i] = secondRecordsFile;

	  ++idx1; ++idx2;

	}
	else if(bamrec_lt(seqs1.recs[idx1], seqs2.recs[idx2])) {
	  seqs1Files[idx1] = first_out;	  
	  ++idx1;
	}
	else {
	  seqs2Files[idx2] = second_out;
	  ++idx2;
	}
		   
      }

      for(;idx1 < seqs1.recs.size(); ++idx1) 
	seqs1Files[idx1] = first_out;	  	

      for(;idx2 < seqs2.recs.size(); ++idx2) 
	seqs2Files[idx2] = second_out;

      // Figure out whether we're splitting the mates up in either case.
      // If they are split up, clear mate information to make the file consistent.

      if(!uniqueValue(seqs1Files))
	clearMateInfo(seqs1);

      if(!uniqueValue(seqs2Files))
	clearMateInfo(seqs2);

      for(int i = 0, ilim = seqs1.recs.size(); i != ilim; ++i) {

	if(seqs1Files[i])
	  seqs1Files[i]->write1(1, seqs1.recs[i]);

      }

      for(int i = 0, ilim = seqs2.recs.size(); i != ilim; ++i) {

	if(seqs2Files[i])
	  seqs2Files[i]->write1(2, seqs2.recs[i]);

      }

    }
    else if(qname_cmp(qname1.c_str(), qname2.c_str()) < 0) {

      if(first_out)
	first_out->write1(1, in1.rec);
      in1.next();

    }
    else {

      if(second_out)
	second_out->write1(2, in2.rec);
      in2.next();

    }

  }

  // One or other file has reached EOF. Write the remainder as first- or second-only records.

  if(first_out) {
    while(!in1.is_eof()) {
      first_out->write1(1, in1.rec);
      in1.next();
    }
  }

  if(second_out) {
    while(!in2.is_eof()) {
      second_out->write1(2, in2.rec);
      in2.next();
    }
  }

  hts_close(in1hf);
  hts_close(in2hf);
  if(first_out)
    htswrapper_close(first_out);
  if(second_out)
    htswrapper_close(second_out);
  if(firstbetter_out)
    htswrapper_close(firstbetter_out);
  if(secondbetter_out)
    htswrapper_close(secondbetter_out);
  if(firstworse_out)
    htswrapper_close(firstworse_out);
  if(secondworse_out)
    htswrapper_close(secondworse_out);

}
